#include <CanvasPoint.h>
#include <CanvasTriangle.h>
#include <Colour.h>
#include <DrawingWindow.h>
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
    float betweenValue = (to - from) / float (numberOfValues - 1);
    for (size_t i = 0; i < numberOfValues; i++) {
        result.push_back(from + i * betweenValue);
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
        result.push_back(from + glm::vec3(i) * betweenValue);
    }
    return result;
}

uint32_t colourPalette(Colour colour) {
    return (255 << 24) + (int(colour.red) << 16) + (int(colour.green) << 8) + int(colour.blue);
}

Colour randomColour() {
    return Colour(rand() % 256, rand() % 256, rand() % 256);
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

void draw(DrawingWindow &window) {

}

/*
// Grayscale
void draw(DrawingWindow &window) {
	window.clearPixels();
	std::vector<float> grayscale = interpolateSingleFloats(255, 0, WIDTH);
	for (size_t y = 0; y < window.height; y++) {
		for (size_t x = 0; x < window.width; x++) {
			float red = grayscale[x];
			float green = grayscale[x];
			float blue = grayscale[x];
			uint32_t colour = (255 << 24) + (int(red) << 16) + (int(green) << 8) + int(blue);
			window.setPixelColour(x, y, colour);
		}
	}
}
*/

/*
Week2
void draw(DrawingWindow &window) {
	window.clearPixels();

	// Task 5: Two Dimensional Colour Interpolation
    glm::vec3 topLeft(255, 0, 0);        // red
    glm::vec3 topRight(0, 0, 255);       // blue
    glm::vec3 bottomRight(0, 255, 0);    // green
    glm::vec3 bottomLeft(255, 255, 0);   // yellow

    // Interpolate colours for the left and right columns
    std::vector<glm::vec3> leftMostColumn = interpolateThreeElementValues(topLeft, bottomLeft, window.height);
    std::vector<glm::vec3> rightMostColumn = interpolateThreeElementValues(topRight, bottomRight, window.height);

	for (size_t y = 0; y < window.height; y++) {
	    // Interpolate colours for each current row
	    std::vector<glm::vec3> eachRowColours = interpolateThreeElementValues(leftMostColumn[y], rightMostColumn[y], window.width);
		for (size_t x = 0; x < window.width; x++) {
		    // Convert glm::vec3 colour to uint32_t
		    uint32_t pixelColour = colourPalette(glm::vec3(eachRowColours[x]));
			// Set the pixel colour using the converted colour
			window.setPixelColour(x, y, pixelColour);

            // Task 3: Single Dimension Greyscale Interpolation
            // Use the x-coordinate as the intensity for a left-to-right gradient
            float grayScale = static_cast<float>(x) / window.width * 255.0;
            // Invert a gradient from white to black
            grayScale = 255.0 - grayScale;
        	// Pack RBG channels into a 32-bit integer
        	uint32_t colour = (255 << 24) + (int(grayScale) << 16) + (int(grayScale) << 8) + int(grayScale);
        	window.setPixelColour(x, y, colour);
		}
	}
}
*/

/*
Week3 Task2
void draw(DrawingWindow &window) {
	window.clearPixels();
	Colour white = Colour(255, 255, 255);
	drawLine(CanvasPoint(0, 0), CanvasPoint(window.width / 2, window.height / 2), window, white);
    drawLine(CanvasPoint(window.width, 0), CanvasPoint(window.width / 2, window.height / 2), window, white);
    drawLine(CanvasPoint(window.width / 2, 0), CanvasPoint(window.width / 2, window.height), window, white);
    drawLine(CanvasPoint(window.width / 3, window.height / 2), CanvasPoint(2 * window.width / 3, window.height / 2), window, white);
}
*/

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
    size_t numberOfVectors = 4;
    std::vector<glm::vec3> resultTEV = interpolateThreeElementValues(from, to, numberOfVectors);
    for (size_t i = 0; i < resultTEV.size(); i++) {
        std::cout << "(" << resultTEV[i].x << ", " << resultTEV[i].y << ", " << resultTEV[i].z << ")" << std::endl;
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
