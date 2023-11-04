#include <CanvasPoint.h>
#include <CanvasTriangle.h>
#include <Colour.h>
#include <DrawingWindow.h>
#include <ModelTriangle.h>
#include <RayTriangleIntersection.h>
#include <TextureMap.h>
#include <TexturePoint.h>
#include <Utils.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <glm/glm.hpp>
#include <map>

#define WIDTH 320
#define HEIGHT 240


bool orbitAnimation = true;
float focalLength = 2.0;
glm::vec3 cameraPosition = glm::vec3(0.0, 0.0, 5.0);
glm::mat3 cameraOrientation = glm::mat3(glm::vec3(1, 0, 0),
                                        glm::vec3(0, 1, 0),
                                        glm::vec3(0, 0, 1));

std::vector<std::vector<float>> initialiseDepthBuffer(int width, int height) {
    std::vector<std::vector<float>> depthBuffer;
    for (int y = 0; y < height; y++) {
        std::vector<float> row;
        for (int x = 0; x < width; x++) {
            row.push_back(INFINITY);
        }
        depthBuffer.push_back(row);
    }
    return depthBuffer;
}

std::vector<float> interpolateSingleFloats (float from, float to, size_t numberOfValues) {
    std::vector<float> result;
    if (numberOfValues == 0 ) {
        result.push_back(from);
        return result;
    }
    float betweenValue = (to - from) / float(numberOfValues);
    for (size_t i = 0; i < numberOfValues; i++) {
        result.push_back(from + float(i) * betweenValue);
    }
    return result;
}

std::vector<glm::vec3> interpolateThreeElementValues(glm::vec3 from, glm::vec3 to, size_t numberOfVectors) {
    std::vector<glm::vec3> result;
    if (numberOfVectors == 0) {
        result.push_back(from);
        return result;
    }
    glm::vec3 betweenValue = (to - from) / glm::vec3(numberOfVectors);
    for (size_t i = 0; i < numberOfVectors; i++) {
        result.push_back(from + betweenValue * glm::vec3(i));
    }
    return result;
}

float interpolation(float targetPosition, float startPosition, float endPosition, float startValue, float endValue) {
    // shouldn't be equal to avoid division by zero
    if (startPosition == endPosition) return startValue;
    // this y = y1 + (x - x1) * ((y2 - y1) / (x2 - x1));
    float slope = (endValue - startValue) / (endPosition - startPosition);
    return startValue + (targetPosition - startPosition) * slope;
}

uint32_t colourPalette(Colour colour) {
    return (255 << 24) + (int(colour.red) << 16) + (int(colour.green) << 8) + int(colour.blue);
}

Colour randomColour() {
    return Colour(rand() % 256, rand() % 256, rand() % 256);
}

CanvasTriangle randomLines() {
    return CanvasTriangle (
            CanvasPoint(rand() % WIDTH, rand() % HEIGHT),
            CanvasPoint(rand() % WIDTH, rand() % HEIGHT),
            CanvasPoint(rand() % WIDTH, rand() % HEIGHT)
            );
}

void drawLine(CanvasPoint from, CanvasPoint to, Colour colour, std::vector<std::vector<float>> &depthBuffer, DrawingWindow &window) {
    uint32_t colourSet = colourPalette(colour);

    float xDiff = to.x - from.x;
    float yDiff = to.y - from.y;
    float zDiff = to.depth - from.depth;

    float numberOfSteps = std::max(abs(xDiff), std::abs(yDiff));

    float xStepSize = xDiff / numberOfSteps;
    float yStepSize = yDiff / numberOfSteps;
    float zStepSize = zDiff / numberOfSteps;

    for (float i= 0.0; i <= numberOfSteps; i++) {
        float x = from.x + (xStepSize * i);
        float y = from.y + (yStepSize * i);
        float z = 1 / (from.depth + (zStepSize * i));

        int roundX = round(x);
        int roundY = round(y);

        if (roundX >= 0 && roundX < WIDTH && roundY >= 0 && roundY < HEIGHT) {
            if (z < depthBuffer[roundY][roundX]) {
                window.setPixelColour(roundX, roundY, colourSet);
                depthBuffer[roundY][roundX] = z;
            }
        }
    }
}

void drawTriangle(CanvasTriangle triangle, Colour colour, std::vector<std::vector<float>> &depthBuffer, DrawingWindow &window) {
    drawLine(triangle.v0(), triangle.v1(), colour, depthBuffer, window);
    drawLine(triangle.v1(), triangle.v2(), colour, depthBuffer, window);
    drawLine(triangle.v2(), triangle.v0(), colour, depthBuffer, window);
}

void drawRandomTriangle(std::vector<std::vector<float>> &depthBuffer, DrawingWindow &window) {
    drawTriangle(randomLines(), randomColour(), depthBuffer, window) ;
}

void triangleRasteriser(CanvasTriangle triangle, Colour colour, std::vector<std::vector<float>> &depthBuffer, DrawingWindow &window) {
    // Sort the vertices by y coordinates
    if (triangle.v0().y > triangle.v1().y) std::swap(triangle.v0(), triangle.v1());
    if (triangle.v0().y > triangle.v2().y) std::swap(triangle.v0(), triangle.v2());
    if (triangle.v1().y > triangle.v2().y) std::swap(triangle.v1(), triangle.v2());

    // Interpolate values for top and bottom part
    std::vector<float> xTop = interpolateSingleFloats(triangle.v0().x, triangle.v2().x, triangle.v1().y - triangle.v0().y);
    std::vector<float> zTop = interpolateSingleFloats(triangle.v0().depth, triangle.v2().depth, triangle.v1().y - triangle.v0().y);
    std::vector<float> xBottom = interpolateSingleFloats(triangle.v1().x, triangle.v2().x, triangle.v2().y - triangle.v1().y);
    std::vector<float> zBottom = interpolateSingleFloats(triangle.v1().depth,  triangle.v2().depth, triangle.v2().y - triangle.v1().y);

    // Extract the x, y coordinates and z coordinates
    float x0 = triangle.v0().x;
    float x1 = triangle.v1().x;
    float x2 = triangle.v2().x;
    float y0 = triangle.v0().y;
    float y1 = triangle.v1().y;
    float y2 = triangle.v2().y;
    float z0 = triangle.v0().depth;
    float z1 = triangle.v1().depth;
    float z2 = triangle.v2().depth;

    // the top part of the triangle
    for (size_t i = 0; i < xTop.size(); i++) {
        float y = triangle.v0().y + i;
        float xStart = xTop[i];
        float xEnd = interpolation(y, y0, y1, x0, x1);
        float zStart = zTop[i];
        float zEnd = interpolation(y, y0, y1, z0, z1);
        drawLine(CanvasPoint(xStart, y, zStart), CanvasPoint(xEnd, y, zEnd), colour,depthBuffer, window);
    }

    // the bottom part of the triangle
    for (size_t i = 0; i < xBottom.size(); i++) {
        float y = triangle.v1().y + i;
        float xStart = xBottom[i];
        float xEnd = interpolation(y, y0, y2, x0, x2);
        float zStart = zBottom[i];
        float zEnd = interpolation(y, y0, y2, z0, z2);
        drawLine(CanvasPoint(xStart, y, zStart), CanvasPoint(xEnd, y, zEnd), colour, depthBuffer, window);
    }
    drawTriangle(triangle, colour, depthBuffer, window);
}

CanvasPoint getCanvasIntersectionPoint(glm::vec3 vertexPosition, float frame) {
    glm::vec3 transposition = (vertexPosition - cameraPosition) * cameraOrientation;
    float u = (-1) * (focalLength * ((transposition.x) / (transposition.z)));
    float v = (focalLength * ((transposition.y) / (transposition.z)));
    // rounding to ensure whole pixels
    CanvasPoint result = CanvasPoint(std::round(u * frame + WIDTH/2.0f), std::round(v * frame + HEIGHT/2.0f));
    result.depth = transposition.z;
    return result;
}

void wireframeRender(std::vector<ModelTriangle> &modelTriangle, std::vector<std::vector<float>> &depthBuffer, DrawingWindow &window) {
    for(auto &triangle : modelTriangle) {
        CanvasTriangle canvasTriangle;
        for(int i = 0; i < 3; i++) {
            glm::vec3 vertex = triangle.vertices[i];
            canvasTriangle.vertices[i] = getCanvasIntersectionPoint(vertex, 240);
        }
        Colour whiteLine = Colour(255, 255, 255);
        drawTriangle(canvasTriangle, whiteLine, depthBuffer, window);
    }
}

void rasterisedRender(std::vector<ModelTriangle> &modelTriangle, std::vector<std::vector<float>> &depthBuffer, DrawingWindow &window) {
    for(auto &triangle: modelTriangle) {
        CanvasTriangle canvasTriangle;
        for(int i = 0; i < 3; i++) {
            glm::vec3 vertex = triangle.vertices[i];
            canvasTriangle.vertices[i] = getCanvasIntersectionPoint(vertex, 240);
        }
        triangleRasteriser(canvasTriangle, triangle.colour, depthBuffer, window);
    }
}

std::map<std::string, Colour> readMTL(const std::string &fileName) {
    std::map<std::string, Colour> readPalette;
    std::ifstream file(fileName);
    std::string line, colourName;
    while (std::getline(file, line)) {
        auto tokens = split(line, ' ');
        if (tokens.empty()) continue;
        if (tokens[0] == "newmtl" && tokens.size() >= 2) {
            colourName = tokens[1];
        } else if (tokens[0] == "Kd" && tokens.size() >= 4) {
            readPalette[colourName] = Colour(colourName, int(std::stof(tokens[1]) * 255), int(std::stof(tokens[2]) * 255), int(std::stof(tokens[3]) * 255));
        }
    }
    file.close();
    return readPalette;
}

std::vector<ModelTriangle> readOBJ(const std::string &fileName, const std::map<std::string, Colour> &readPalette, float scalingFactor) {
    std::vector<ModelTriangle> triangles;
    std::vector<glm::vec3> vertices;
    std::ifstream file(fileName);
    std::string line;
    Colour colour;

    while(std::getline(file, line)) {
        std::vector<std::string> tokens = split(line, ' ');
        if (tokens.empty()) continue;
        if (tokens[0] == "v" && tokens.size() >= 4) {
            vertices.push_back(glm::vec3 (
                    std::stof(tokens[1]) * scalingFactor,
                    std::stof(tokens[2]) * scalingFactor,
                    std::stof(tokens[3]) * scalingFactor));
        } else if (tokens[0] == "usemtl" && tokens.size() >= 2) {
            colour = readPalette.at(tokens[1]);
        } else if (tokens[0] == "f" && tokens.size() >= 4) {
            triangles.push_back(ModelTriangle(vertices[std::stoi(tokens[1]) - 1],
                                              vertices[std::stoi(tokens[2]) - 1],
                                              vertices[std::stoi(tokens[3]) - 1],
                                              colour));
        }
    }
    file.close();
    return triangles;
}

void lookAt() {
    glm::vec3 forward = glm::normalize(cameraPosition);
    glm::vec3 right = glm::normalize(glm::cross(glm::vec3(0, 1, 0), forward));
    glm::vec3 up = glm::cross(forward, right);
    cameraOrientation = glm::mat3(right, up, forward);
}

void orbitCamera() {
    float speed = 1.0f;
    float theta = glm::radians(speed);
    glm::mat3 rotationMatrix = glm::mat3(glm::vec3(glm::cos(theta), 0, glm::sin(theta)),
                                    glm::vec3(0, 1, 0),
                                    glm::vec3(-glm::sin(theta), 0, glm::cos(theta)));
    cameraPosition = rotationMatrix * cameraPosition;
    lookAt();
}

void draw(DrawingWindow &window) {
    window.clearPixels();
    std::vector<std::vector<float>> depthBuffer = initialiseDepthBuffer(WIDTH, HEIGHT);
    std::map<std::string, Colour> readPalette = readMTL("./src/cornell-box.mtl");
    std::vector<ModelTriangle> modelTriangle = readOBJ("./src/cornell-box.obj", readPalette, 0.35);
    rasterisedRender(modelTriangle, depthBuffer, window);
    if (orbitAnimation) {
        orbitCamera();
    }
    // rasterisedRender(modelTriangle, depthBuffer, window);
}

void handleEvent(SDL_Event event, std::vector<std::vector<float>> &depthBuffer, DrawingWindow &window) {
    if (event.type == SDL_KEYDOWN) {
        bool cameraUpdated = false;
        float move = 0.1f;
        float theta = M_PI / 2.0f;
        if (event.key.keysym.sym == SDLK_LEFT) {
            cameraPosition.x += move;
            std::cout << "LEFT" << std::endl;
            cameraUpdated = true;
        }
        else if (event.key.keysym.sym == SDLK_RIGHT) {
            cameraPosition.x -= move;
            std::cout << "RIGHT" << std::endl;
            cameraUpdated = true;
        }
        else if (event.key.keysym.sym == SDLK_UP) {
            cameraPosition.y += move;
            std::cout << "UP" << std::endl;
            cameraUpdated = true;
        }
        else if (event.key.keysym.sym == SDLK_DOWN) {
            cameraPosition.y -= move;
            std::cout << "DOWN" << std::endl;
            cameraUpdated = true;
        }
        else if (event.key.keysym.sym == SDLK_i) {
            cameraPosition = glm::mat3(glm::vec3(1, 0, 0),
                                          glm::vec3(0, glm::cos(theta), -glm::sin(theta)),
                                          glm::vec3(0, glm::sin(theta), glm::cos(theta))) * cameraPosition;
            cameraUpdated = true;
        }
        else if (event.key.keysym.sym == SDLK_j) {
            cameraPosition = glm::mat3(glm::vec3(glm::cos(theta), 0, glm::sin(theta)),
                                          glm::vec3(0, 1, 0),
                                          glm::vec3(-glm::sin(theta), 0, glm::cos(theta))) * cameraPosition;
            cameraUpdated = true;
        }
        else if (event.key.keysym.sym == SDLK_k) {
            cameraPosition = glm::mat3(glm::vec3(1, 0, 0),
                                          glm::vec3(0, glm::cos(-theta), -glm::sin(-theta)),
                                          glm::vec3(0, glm::sin(-theta), glm::cos(-theta))) * cameraPosition;
            cameraUpdated = true;
        }
        else if (event.key.keysym.sym == SDLK_l) {
            cameraPosition = glm::mat3(glm::vec3(glm::cos(-theta), 0, glm::sin(-theta)),
                                          glm::vec3(0, 1, 0),
                                          glm::vec3(-glm::sin(-theta), 0, glm::cos(-theta))) * cameraPosition;
            cameraUpdated = true;
        }
        else if (event.key.keysym.sym == SDLK_o) {
            orbitAnimation = !orbitAnimation;
            std::cout << "Orbit Camera " << (orbitAnimation ? "ON" : "OFF") << std::endl;
        }
        if (cameraUpdated) {
            lookAt();
            // draw(window);
        }
    } else if (event.type == SDL_MOUSEBUTTONDOWN) {
        window.savePPM("output.ppm");
        window.saveBMP("output.bmp");
    }
}

int main(int argc, char *argv[]) {
    DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);

    std::vector<std::vector<float>> depthBuffer = initialiseDepthBuffer(WIDTH, HEIGHT);

    SDL_Event event;

    while (true) {
        // We MUST poll for events - otherwise the window will freeze !
        if (window.pollForInputEvents(event)) handleEvent(event, depthBuffer, window);
        draw(window);
        // Need to render the frame at the end, or nothing actually gets shown on the screen !
        window.renderFrame();
    }
}
