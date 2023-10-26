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
float focalLength = 2.0;

std::vector<float> interpolateSingleFloats (float from, float to, size_t numberOfValues) {
    std::vector<float> result;
    if (numberOfValues == 0 ) {
        result.push_back(from);
        return result;
    }
    float betweenValue = (to - from) / float(numberOfValues - 1);
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
    glm::vec3 betweenValue = (to - from) / glm::vec3(numberOfVectors - 1);
    for (size_t i = 0; i < numberOfVectors; i++) {
        result.push_back(from + betweenValue * glm::vec3(i));
    }
    return result;
}

uint32_t colourPalette(Colour colour) {
    return (255 << 24) + (int(colour.red) << 16) + (int(colour.green) << 8) + int(colour.blue);
}

Colour randomColour() {
    return Colour(rand() % 256, rand() % 256, rand() % 256);
}

float interpolation(float targetPosition, float startPosition, float endPosition, float startValue, float endValue) {
    // shouldn't be equal to avoid division by zero
    if (startPosition == endPosition) return startValue;

    // this y = y1 + (x - x1) * ((y2 - y1) / (x2 - x1));
    float slope = (endValue - startValue) / (endPosition - startPosition);
    float targetValue = startValue + (targetPosition - startPosition) * slope;
    return targetValue;
}

CanvasTriangle randomLines() {
    return CanvasTriangle (
            CanvasPoint(rand() % WIDTH, rand() % HEIGHT),
            CanvasPoint(rand() % WIDTH, rand() % HEIGHT),
            CanvasPoint(rand() % WIDTH, rand() % HEIGHT)
            );
}

void drawLine(CanvasPoint from, CanvasPoint to, Colour colour, DrawingWindow &window) {
    uint32_t colourSet = colourPalette(colour);
    float xDiff = to.x - from.x;
    float yDiff = to.y - from.y;
    float numberOfSteps = std::max(abs(xDiff), std::abs(yDiff));

    float xStepSize = xDiff / numberOfSteps;
    float yStepSize = yDiff / numberOfSteps;

    for (float i= 0.0; i <= numberOfSteps; i++) {
        float x = from.x + (xStepSize * i);
        float y = from.y + (yStepSize * i);

        if (x >= 0 && x < WIDTH && y >= 0 && y < HEIGHT) {
                window.setPixelColour(round(x), round(y), colourSet);
        }
    }
}

void drawTriangle(CanvasTriangle triangle, Colour colour, DrawingWindow &window) {
    drawLine(triangle.v0(), triangle.v1(), colour, window);
    drawLine(triangle.v1(), triangle.v2(), colour, window);
    drawLine(triangle.v2(), triangle.v0(), colour, window);
}

void drawRandomTriangle(DrawingWindow &window) {
    drawTriangle(randomLines(), randomColour(), window) ;
}

void triangleRasteriser(CanvasTriangle triangle, Colour colour, DrawingWindow &window) {
    // Sort the vertices by y coordinates
    if (triangle.v0().y > triangle.v1().y) std::swap(triangle.v0(), triangle.v1());
    if (triangle.v0().y > triangle.v2().y) std::swap(triangle.v0(), triangle.v2());
    if (triangle.v1().y > triangle.v2().y) std::swap(triangle.v1(), triangle.v2());

    // Extract the x and y coordinates
    float x0 = triangle.v0().x;
    float x1 = triangle.v1().x;
    float x2 = triangle.v2().x;
    float y0 = triangle.v0().y;
    float y1 = triangle.v1().y;
    float y2 = triangle.v2().y;


    // Calculate slopes
    float slope01 = (y1 - y0) != 0 ? (x1 - x0) / (y1 - y0) : 0;
    float slope02 = (y2 - y0) != 0 ? (x2 - x0) / (y2 - y0) : 0;
    float slope12 = (y2 - y1) != 0 ? (x2 - x1) / (y2 - y1) : 0;


    // the top part of the triangle
    for (size_t y = y0; y <= y1; y++) {
        float xStart = x0 + (y - y0) * slope02;
        float xEnd = x0 + (y - y0) * slope01;
        drawLine(CanvasPoint(xStart, y), CanvasPoint(xEnd, y), colour,window);
    }

    // the bottom part of the triangle
    for (size_t y = y1 + 1; y <= y2; y++) {
        float xStart = x1 + (y - y1) * slope12;
        float xEnd = x0 + (y - y0) * slope02;
        drawLine(CanvasPoint(xStart, y), CanvasPoint(xEnd, y), colour, window);
    }
    drawTriangle(triangle, colour, window);
}

CanvasPoint getCanvasIntersectionPoint(glm::vec3 vertexPosition, float frame) {
    glm::vec3 transposition = vertexPosition - cameraPosition;
    float u = (-1) * (focalLength * ((transposition.x) / (transposition.z)));
    float v = (focalLength * ((transposition.y) / (transposition.z)));
    // rounding to ensure whole pixels
    CanvasPoint result = CanvasPoint(round(u * frame + WIDTH/2.0f), round(v * frame + HEIGHT/2.0f));
    return result;
}

void wireframeRender(std::vector<ModelTriangle> &modelTriangle, DrawingWindow &window) {
    for(auto &triangle : modelTriangle) {
        CanvasTriangle canvasTriangle;
        for(int i = 0; i < 3; i++) {
            glm::vec3 vertex = triangle.vertices[i];
            canvasTriangle.vertices[i] = getCanvasIntersectionPoint(vertex, 240);
        }
        Colour whiteLine = Colour(255, 255, 255);
        drawTriangle(canvasTriangle, whiteLine, window);
    }
}

void rasterisedRender(std::vector<ModelTriangle> &modelTriangle, DrawingWindow &window) {
    for(auto &triangle: modelTriangle) {
        CanvasTriangle canvasTriangle;
        for(int i = 0; i < 3; i++) {
            glm::vec3 vertex = triangle.vertices[i];
            canvasTriangle.vertices[i] = getCanvasIntersectionPoint(vertex, 240);
        }
        triangleRasteriser(canvasTriangle, triangle.colour, window);
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

void handleEvent(SDL_Event event, DrawingWindow &window) {
    if (event.type == SDL_KEYDOWN) {
        if (event.key.keysym.sym == SDLK_LEFT) std::cout << "LEFT" << std::endl;
        else if (event.key.keysym.sym == SDLK_RIGHT) std::cout << "RIGHT" << std::endl;
        else if (event.key.keysym.sym == SDLK_UP) std::cout << "UP" << std::endl;
        else if (event.key.keysym.sym == SDLK_DOWN) std::cout << "DOWN" << std::endl;
        else if (event.key.keysym.sym == SDLK_u) {
            drawRandomTriangle(window);
            window.renderFrame();
        }
        else if (event.key.keysym.sym == SDLK_f) {
            triangleRasteriser(randomLines(), randomColour(), window);
            window.renderFrame();
        }
    } else if (event.type == SDL_MOUSEBUTTONDOWN) {
        window.savePPM("output.ppm");
        window.saveBMP("output.bmp");
    }
}

int main(int argc, char *argv[]) {
    DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);

    std::map<std::string, Colour> readPalette = readMTL("./src/cornell-box.mtl");
    std::vector<ModelTriangle> modelTriangle = readOBJ("./src/cornell-box.obj", readPalette, 0.35);

    // wireframeRender(modelTriangle, window);
    rasterisedRender(modelTriangle, window);

    SDL_Event event;

    while (true) {
        // We MUST poll for events - otherwise the window will freeze !
        if (window.pollForInputEvents(event)) handleEvent(event, window);
        // Need to render the frame at the end, or nothing actually gets shown on the screen !
        window.renderFrame();
    }
}
