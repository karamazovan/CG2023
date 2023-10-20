#include <CanvasPoint.h>
#include <CanvasTriangle.h>
#include <Colour.h>
#include <DrawingWindow.h>
#include <TextureMap.h>
#include <TexturePoint.h>
#include <Utils.h>
#include <fstream>
#include <iostream>
#include <vector>
#include <glm/glm.hpp>


#define WIDTH 320
#define HEIGHT 240

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
    glm::vec3 divisor(numberOfVectors - 1, numberOfVectors - 1, numberOfVectors - 1);
    glm::vec3 betweenValue = (to - from) / divisor;
    // glm::vec3 betweenValue = (to - from) / glm::vec3(numberOfVectors - 1);
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

float interpolation(float x, float x1, float x2, float y1, float y2) {
    // shouldn't be equal to avoid division by zero
    if (x1 == x2) {
        return y1;
    }
    float slope = (y2 - y1) / (x2 - x1);
    float y = y1 + (x - x1) * slope;
    return y;
}

CanvasTriangle randomLines() {
    CanvasPoint v0 = CanvasPoint(rand() % WIDTH, rand() % HEIGHT);
    CanvasPoint v1 = CanvasPoint(rand() % WIDTH, rand() % HEIGHT);
    CanvasPoint v2 = CanvasPoint(rand() % WIDTH, rand() % HEIGHT);
    return CanvasTriangle(v0, v1, v2);
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

void drawLineTexture(CanvasPoint from, CanvasPoint to, TextureMap &textureMap, DrawingWindow &window) {
    float xDiff = to.x - from.x;
    float yDiff = to.y - from.y;
    float numberOfSteps = std::max(abs(xDiff), std::abs(yDiff));

    float xStepSize = xDiff / numberOfSteps;
    float yStepSize = yDiff / numberOfSteps;

    float xTextureDiff = to.texturePoint.x - from.texturePoint.x;
    float yTextureDiff = to.texturePoint.y - from.texturePoint.y;
    float xTextureStepSize = xTextureDiff / numberOfSteps;
    float yTextureStepSize = yTextureDiff / numberOfSteps;

    for (float i = 0.0; i < numberOfSteps; i++) {
        float x = std::round(from.x + (xStepSize * i));
        float y = std::round(from.y + (yStepSize * i));
        float xTexture = from.texturePoint.x + xTextureStepSize * i;
        float yTexture = from.texturePoint.y + yTextureStepSize * i;

        if (x >= 0 && x < WIDTH && y >= 0 && y < HEIGHT && xTexture >= 0 && xTexture < textureMap.width && yTexture >= 0 && yTexture < textureMap.height) {
            uint32_t textureColour = textureMap.pixels[size_t(xTexture) + size_t(yTexture) * textureMap.width];
            window.setPixelColour(round(x), round(y), textureColour);
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
        drawLine(CanvasPoint(xStart, y), CanvasPoint(xEnd, y), colour, window);
    }

    // the bottom part of the triangle
    for (size_t y = y1 + 1; y <= y2; y++) {
        float xStart = x1 + (y - y1) * slope12;
        float xEnd = x0 + (y - y0) * slope02;
        drawLine(CanvasPoint(xStart, y), CanvasPoint(xEnd, y), colour, window);
    }
    // draw white line
    Colour whiteLine = Colour(255, 255, 255);
    drawTriangle(triangle, whiteLine, window);
}

void textureMapper(CanvasTriangle triangle, TextureMap &textureMap, DrawingWindow &window) {
    // Sort the vertices by y coordinates
    if (triangle.v0().y > triangle.v1().y) std::swap(triangle.v0(), triangle.v1());
    if (triangle.v0().y > triangle.v2().y) std::swap(triangle.v0(), triangle.v2());
    if (triangle.v1().y > triangle.v2().y) std::swap(triangle.v1(), triangle.v2());

    for (size_t y = triangle.v0().y; y < triangle.v1().y; y++) {
        CanvasPoint leftPoint, rightPoint;

        leftPoint.x = interpolation(y, triangle.v0().y, triangle.v1().y, triangle.v0().x, triangle.v1().x);
        leftPoint.texturePoint.x = interpolation(y,triangle.v0().y, triangle.v1().y, triangle.v0().texturePoint.x, triangle.v1().texturePoint.x);
        leftPoint.texturePoint.y = interpolation(y,triangle.v0().y, triangle.v1().y, triangle.v0().texturePoint.y, triangle.v1().texturePoint.y);

        rightPoint.x = interpolation(y, triangle.v0().y, triangle.v2().y, triangle.v0().x, triangle.v2().x);
        rightPoint.texturePoint.x = interpolation(y,triangle.v0().y, triangle.v2().y, triangle.v0().texturePoint.x, triangle.v2().texturePoint.x);
        rightPoint.texturePoint.y = interpolation(y,triangle.v0().y, triangle.v2().y, triangle.v0().texturePoint.y, triangle.v2().texturePoint.y);

        leftPoint.y = rightPoint.y = y;

        drawLineTexture(leftPoint, rightPoint, textureMap, window);
    }

    for (size_t y = triangle.v1().y; y < triangle.v2().y; y++) {
        CanvasPoint leftPoint, rightPoint;

        leftPoint.x = interpolation(y, triangle.v1().y, triangle.v2().y, triangle.v1().x, triangle.v2().x);
        leftPoint.texturePoint.x = interpolation(y,triangle.v1().y, triangle.v2().y, triangle.v1().texturePoint.x, triangle.v2().texturePoint.x);
        leftPoint.texturePoint.y = interpolation(y,triangle.v1().y, triangle.v2().y, triangle.v1().texturePoint.y, triangle.v2().texturePoint.y);

        rightPoint.x = interpolation(y, triangle.v0().y, triangle.v2().y, triangle.v0().x, triangle.v2().x);
        rightPoint.texturePoint.x = interpolation(y,triangle.v0().y, triangle.v2().y, triangle.v0().texturePoint.x, triangle.v2().texturePoint.x);
        rightPoint.texturePoint.y = interpolation(y,triangle.v0().y, triangle.v2().y, triangle.v0().texturePoint.y, triangle.v2().texturePoint.y);

        leftPoint.y = rightPoint.y = y;

        drawLineTexture(leftPoint, rightPoint, textureMap, window);
    }
    Colour whiteLine = Colour(255, 255, 255);
    drawTriangle(triangle, whiteLine, window);
}

void draw(DrawingWindow &window) {
    window.clearPixels();
    TextureMap textureMap("./src/texture.ppm");
    CanvasTriangle triangle(CanvasPoint(160, 10), CanvasPoint(300, 230), CanvasPoint(10, 150));
    triangle.v0().texturePoint = TexturePoint(195, 5);
    triangle.v1().texturePoint = TexturePoint(395, 380);
    triangle.v2().texturePoint = TexturePoint(65, 330);
    textureMapper(triangle, textureMap, window);
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

    // W2 - Task 2: Single Element Numerical Interpolation
	std::vector<float> result;
    result = interpolateSingleFloats(2.2, 8.5, 7);
	for (size_t i = 0; i < result.size(); i++) {
	    std::cout << result[i] << " ";
	}
    std::cout << std::endl;

    // W2 - Task 4: Three Element Numerical Interpolation
    glm::vec3 from(1.0, 4.0, 9.2);
    glm::vec3 to(4.0, 1.0, 9.8);
    for (glm::vec3 vec : interpolateThreeElementValues(from, to, 4)) {
        std::cout << "(" << vec.x << ", " << vec.y << ", " << vec.z << ")" << std::endl;
    }

    SDL_Event event;

    while (true) {
		// We MUST poll for events - otherwise the window will freeze !
		if (window.pollForInputEvents(event)) handleEvent(event, window);
		draw(window);
		// Need to render the frame at the end, or nothing actually gets shown on the screen !
		window.renderFrame();
	}
}
