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

glm::vec3 cameraPosition = glm::vec3(0.0, 0.0, 5.0);
glm::mat3 cameraOrientation = glm::mat3(1, 0, 0, 0, 1, 0, 0, 0, 1);
float focalLength = 2.0;

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

// (x, y, z)
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
    for (size_t y = y0; y <= y1; y++) {
        float xStart = interpolation(y, y0, y2, x0, x2);
        float xEnd = interpolation(y, y0, y1, x0, x1);
        float zStart = interpolation(y, y0, y2, z0, z2);
        float zEnd = interpolation(y, y0, y1, z0, z1);
        drawLine(CanvasPoint(xStart, y, zStart), CanvasPoint(xEnd, y, zEnd), colour,depthBuffer, window);
    }

    // the bottom part of the triangle
    for (size_t y = y1 + 1; y <= y2; y++) {
        float xStart = interpolation(y, y1, y2, x1, x2);
        float xEnd = interpolation(y, y0, y2, x0, x2);
        float zStart = interpolation(y, y1, y2, z1, z2);
        float zEnd = interpolation(y, y0, y2, z0, z2);
        drawLine(CanvasPoint(xStart, y, zStart), CanvasPoint(xEnd, y, zEnd), colour, depthBuffer, window);
    }
    // drawTriangle(triangle, colour, depthBuffer, window);
}

CanvasPoint getCanvasIntersectionPoint(glm::vec3 vertexPosition, float frame) {
    glm::vec3 transposition = vertexPosition - cameraPosition;
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

void rasterisedRender(std::vector<ModelTriangle> &modelTriangle, DrawingWindow &window) {
    window.clearPixels();
    std::vector<std::vector<float>> depthBuffer = initialiseDepthBuffer(WIDTH, HEIGHT);
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

void draw(std::vector<std::vector<float>> &depthBuffer, DrawingWindow &window) {
    window.clearPixels();
}

void handleEvent(SDL_Event event, std::vector<std::vector<float>> &depthBuffer, DrawingWindow &window) {
    float move = 0.1f;
    float angle = glm::radians(5.0f);
    if (event.type == SDL_KEYDOWN) {
        if (event.key.keysym.sym == SDLK_LEFT) {
            cameraPosition.x += move;
            std::cout << "LEFT" << std::endl;
        }
        else if (event.key.keysym.sym == SDLK_RIGHT) {
            cameraPosition.x -= move;
            std::cout << "RIGHT" << std::endl;
        }
        else if (event.key.keysym.sym == SDLK_UP) {
            cameraPosition.y += move;
            std::cout << "UP" << std::endl;
        }
        else if (event.key.keysym.sym == SDLK_DOWN) {
            cameraPosition.y -= move;
            std::cout << "DOWN" << std::endl;
        }
        else if (event.key.keysym.sym == SDLK_i) {
            cameraPosition = glm::mat3(glm::vec3(1, 0, 0),
                                       glm::vec3(0, cos(-angle), -sin(-angle)),
                                       glm::vec3(0, sin(-angle), cos(-angle))) * cameraPosition;
        }
        else if (event.key.keysym.sym == SDLK_j) {
            cameraPosition = glm::mat3(glm::vec3(cos(-angle), 0, sin(-angle)),
                                       glm::vec3(0, 1, 0),
                                       glm::vec3(-sin(-angle), 0, cos(-angle))) * cameraPosition;
        }
        else if (event.key.keysym.sym == SDLK_k) {
            cameraPosition = glm::mat3(glm::vec3(1, 0, 0),
                                       glm::vec3(0, cos(angle), -sin(angle)),
                                       glm::vec3(0, sin(angle), cos(angle))) * cameraPosition;
        }
        else if (event.key.keysym.sym == SDLK_l) {
            cameraPosition = glm::mat3(glm::vec3(cos(angle), 0, sin(angle)),
                                       glm::vec3(0, 1, 0),
                                       glm::vec3(-sin(angle), 0, cos(angle))) * cameraPosition;
        }
    } else if (event.type == SDL_MOUSEBUTTONDOWN) {
        window.savePPM("output.ppm");
        window.saveBMP("output.bmp");
    }
}

int main(int argc, char *argv[]) {
    DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);

    std::vector<std::vector<float>> depthBuffer = initialiseDepthBuffer(WIDTH, HEIGHT);

    std::map<std::string, Colour> readPalette = readMTL("./src/cornell-box.mtl");
    std::vector<ModelTriangle> modelTriangle = readOBJ("./src/cornell-box.obj", readPalette, 0.35);


    SDL_Event event;

    while (true) {
        // We MUST poll for events - otherwise the window will freeze !
        if (window.pollForInputEvents(event)) handleEvent(event, depthBuffer, window);
        // draw(depthBuffer, window);
        rasterisedRender(modelTriangle, window);
        // Need to render the frame at the end, or nothing actually gets shown on the screen !
        window.renderFrame();
    }
}
