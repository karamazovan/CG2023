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
#include <algorithm>
#include <unordered_map>

#define WIDTH 320
#define HEIGHT 240

int callDraw = 0;
float focalLength = 2.0f;
bool orbitAnimation = true;
glm::vec3 ambientLighting(0.2f, 0.2f, 0.2f);
glm::mat3 cameraOrientation = glm::mat3(1.0);
glm::vec3 cameraPosition = glm::vec3(0.0f, 0.7f, 2.5f);
glm::vec3 lightPosition = glm::vec3(0.0f, 0.5f, 0.5f);

int roundToInt(float value) {
    return int(std::round(value));
}

glm::vec3 modelTriangleNormal(ModelTriangle &triangle) {
    return glm::normalize(glm::cross(triangle.vertices[1] - triangle.vertices[0], triangle.vertices[2] - triangle.vertices[0]));
}

std::vector<std::vector<float>> initialiseDepthBuffer(int width, int height) {
    std::vector<std::vector<float>> depthBuffer;
    for (int y = 0; y < height; y++) {
        std::vector<float> row;
        for (int x = 0; x < width; x++) {
            row.push_back(0);
        }
        depthBuffer.push_back(row);
    }
    return depthBuffer;
}

std::vector<float> interpolateSingleFloats(float from, float to, size_t numberOfValues) {
    std::vector<float> result;
    if (numberOfValues == 0) {
        result.push_back(from);
        return result;
    }
    float betweenValue = (to - from) / float(numberOfValues - 1);
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
    glm::vec3 betweenValue = (to - from) / glm::vec3(numberOfVectors - 1);
    for (size_t i = 0; i < numberOfVectors; i++) {
        result.push_back(from + betweenValue * glm::vec3(i));
    }
    return result;
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

    for (float i = 0.0; i <= numberOfSteps; i++) {
        float x = from.x + (xStepSize * i);
        float y = from.y + (yStepSize * i);
        float z = 1 / (from.depth + (zStepSize * i));

        if (roundToInt(x)>= 0 && roundToInt(x) < WIDTH && roundToInt(y) >= 0 && roundToInt(y)< HEIGHT) {
            if (z > depthBuffer[roundToInt(y)][roundToInt(x)]) {
                window.setPixelColour(roundToInt(x), roundToInt(y), colourSet);
               depthBuffer[roundToInt(y)][roundToInt(x)] = z;
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
    std::vector<CanvasPoint> vertices = {triangle.v0(), triangle.v1(), triangle.v2()};
    std::sort(vertices.begin(), vertices.end(), [](CanvasPoint &a, CanvasPoint &b) {
       return a.y < b.y;
    });

    CanvasPoint &v0 = vertices[0];
    CanvasPoint &v1 = vertices[1];
    CanvasPoint &v2 = vertices[2];

    // Interpolate the x, y and z values
    auto xLeftSide = interpolateSingleFloats(v0.x, v1.x, roundToInt(v1.y) - roundToInt(v0.y) + 1);
    auto zLeftSide = interpolateSingleFloats(v0.depth, v1.depth, roundToInt(v1.y) - roundToInt(v0.y) + 1);
    auto xRightSide = interpolateSingleFloats(v0.x, v2.x, roundToInt(v2.y) - roundToInt(v0.y) + 1);
    auto zRightSide = interpolateSingleFloats(v0.depth, v2.depth, roundToInt(v2.y) - roundToInt(v0.y) + 1);

    // Draw the first half of the triangle
    int yStart = roundToInt(v0.y);
    int yMiddle = roundToInt(v1.y);
    for (int y = yStart; y <= yMiddle; y++) {
        CanvasPoint start(xLeftSide[y - yStart], y, zLeftSide[y - yStart]);
        CanvasPoint end(xRightSide[y - yStart], y, zRightSide[y - yStart]);
        drawLine(start, end, colour, depthBuffer, window);
    }

    // Interpolation for the 2nd half of the triangle
    xLeftSide = interpolateSingleFloats(v1.x, v2.x, roundToInt(v2.y) - roundToInt(v1.y) + 1);
    zLeftSide = interpolateSingleFloats(v1.depth, v2.depth, roundToInt(v2.y) - roundToInt(v1.y) + 1);

    // Draw the 2nd half of the triangle
    int yEnd = roundToInt(v2.y);
    for (int y = yMiddle + 1; y <= yEnd; y++) {
        CanvasPoint start(xLeftSide[y - yMiddle - 1], y, zLeftSide[y - yMiddle - 1]);
        CanvasPoint end(xRightSide[y - yStart], y, zRightSide[y - yStart]);
        drawLine(start, end, colour, depthBuffer, window);
    }
    drawTriangle(triangle, colour, depthBuffer, window);
}

CanvasPoint getCanvasIntersectionPoint(glm::vec3 vertexPosition, float frame) {
    glm::vec3 transposition = (vertexPosition - cameraPosition) * cameraOrientation;
    float u = (-1) * (focalLength * ((transposition.x) / (transposition.z)));
    float v = (focalLength * ((transposition.y) / (transposition.z)));
    // rounding to ensure whole pixels
    CanvasPoint result = CanvasPoint(std::round(u * frame + WIDTH/2.0f), std::round(v * frame + HEIGHT/2.0f));
    result.depth = -transposition.z;
    return result;
}

void wireframeRender(std::vector<ModelTriangle> &modelTriangle, std::vector<std::vector<float>> &depthBuffer, DrawingWindow &window) {
    for (auto &triangle : modelTriangle) {
        CanvasTriangle canvasTriangle;
        for (int i = 0; i < 3; i++) {
            glm::vec3 vertex = triangle.vertices[i];
            canvasTriangle.vertices[i] = getCanvasIntersectionPoint(vertex, 240);
        }
        Colour whiteLine = Colour(255, 255, 255);
        drawTriangle(canvasTriangle, whiteLine, depthBuffer, window);
    }
}

void rasterisedRender(std::vector<ModelTriangle> &modelTriangle, std::vector<std::vector<float>> &depthBuffer, DrawingWindow &window) {
    for (auto &triangle: modelTriangle) {
        CanvasTriangle canvasTriangle;
        for (int i = 0; i < 3; i++) {
            glm::vec3 vertex = triangle.vertices[i];
            canvasTriangle.vertices[i] = getCanvasIntersectionPoint(vertex, 240);
        }
        triangleRasteriser(canvasTriangle, triangle.colour, depthBuffer, window);
    }
}

void lookAt() {
    glm::vec3 forward = glm::normalize(cameraPosition);
    glm::vec3 right = glm::normalize(glm::cross(glm::vec3(0, 1, 0), forward));
    glm::vec3 up = glm::cross(forward, right);
    cameraOrientation = glm::mat3(right, up, forward);
}

void orbitCamera() {
    float speed = 0.5f;
    float theta = glm::radians(speed);
    glm::mat3 rotationMatrix = glm::mat3(glm::vec3(glm::cos(theta), 0, glm::sin(theta)),
                                         glm::vec3(0, 1, 0),
                                         glm::vec3(-glm::sin(theta), 0, glm::cos(theta)));
    cameraPosition = rotationMatrix * cameraPosition;
    lookAt();
}

void cameraTranslation(float tx, float ty, float tz) {
    cameraPosition.x += tx;
    cameraPosition.y += ty;
    cameraPosition.z += tz;
}

void cameraRotation(float angleX, float angleY) {
    if (angleX != 0) {
        float thetaX = glm::radians(angleX);
        glm::mat3 rotationMatrixX = glm::mat3(glm::vec3(1, 0, 0),
                                              glm::vec3(0, glm::cos(thetaX), -glm::sin(thetaX)),
                                              glm::vec3(0, glm::sin(thetaX), glm::cos(thetaX)));
        cameraPosition = rotationMatrixX * cameraPosition;
    }
    if (angleY != 0) {
        float thetaY = glm::radians(angleY);
        glm::mat3 rotationMatrixY = glm::mat3(glm::vec3(glm::cos(thetaY), 0, glm::sin(thetaY)),
                                              glm::vec3( 0, 1, 0),
                                              glm::vec3(-glm::sin(thetaY), 0, glm::cos(thetaY)));
        cameraPosition = rotationMatrixY * cameraPosition;
    }
    lookAt();
}

RayTriangleIntersection getClosestValidIntersection(std::vector<ModelTriangle> &modelTriangle, glm::vec3 &cameraPosition, glm::vec3 rayDirection) {
    RayTriangleIntersection closestIntersection;
    closestIntersection.distanceFromCamera = INFINITY;

    for (size_t i = 0; i < modelTriangle.size(); i++) {
        ModelTriangle triangle = modelTriangle[i];
        glm::vec3 e0 = triangle.vertices[1] - triangle.vertices[0];
        glm::vec3 e1 = triangle.vertices[2] - triangle.vertices[0];
        glm::vec3 SPVector = cameraPosition - triangle.vertices[0];
        glm::mat3 DEMatrix(-rayDirection, e0, e1);
        glm::vec3 possibleSolution = glm::inverse(DEMatrix) * SPVector;

        float t = possibleSolution[0];
        float u = possibleSolution[1];
        float v = possibleSolution[2];

        if (u >= 0.0 && u <= 1.0 && v >= 0.0 && v <= 1.0 && (u + v) <= 1.0) {
            if (t > 0 && t < closestIntersection.distanceFromCamera) {
                closestIntersection.distanceFromCamera = t;
                glm::vec3 r = cameraPosition + t * rayDirection;
                closestIntersection = RayTriangleIntersection(r, t, triangle, i);
            }
        }
    }
    return closestIntersection;
}

float proximityLighting(glm::vec3 point) {
    float lightDistance = glm::length(lightPosition - point);
    float lightIntensity = 10.0f / (4.0f * M_PI * lightDistance * lightDistance);
    return glm::clamp(lightIntensity, 0.0f, 1.0f);
}

float angleOfIncidenceLighting(glm::vec3 point, glm::vec3 normal) {
    glm::vec3 lightDirection = glm::normalize(lightPosition - point);
    float angleIncidence = glm::dot(normal, lightDirection);
    float diffuseIntensity = std::max(angleIncidence, 0.0f);
    return diffuseIntensity;
}

float specularLighting(glm::vec3 point, glm::vec3 normal) {
    glm::vec3 cameraDirection = glm::normalize(cameraPosition - point);
    glm::vec3 lightDirection = glm::normalize(lightPosition - point);
    glm::vec3 reflectionDirection = glm::reflect(-lightDirection, normal);
    float specularIntensity = glm::pow(std::max(glm::dot(cameraDirection, reflectionDirection), 0.0f), 64);
    return specularIntensity * 0.5f;
}

glm::vec3 pixelToDirection(int x, int y) {
    float xNormalise = (x - WIDTH / 2.0f) * (1.0f / (HEIGHT * 2.0f/3.0f));
    float yNormalise = (y - HEIGHT / 2.0f) * (1.0f / (HEIGHT * 2.0f/3.0f));
    glm::vec3 pixelPosition = cameraPosition + cameraOrientation * glm::vec3(xNormalise, -yNormalise, -focalLength);
    return glm::normalize(pixelPosition - cameraPosition);
}

void rasterisedScene(std::vector<ModelTriangle> &modelTriangle,  std::vector<std::vector<float>> &depthBuffer, DrawingWindow &window) {
    for (int y = 0; y < HEIGHT; y++) {
        for (int x = 0; x < WIDTH; x++) {
            // converting from 2D pixel into 3D direction
            glm::vec3 rayDirection = pixelToDirection(x, y);
            RayTriangleIntersection intersectioPoint = getClosestValidIntersection(modelTriangle, cameraPosition, rayDirection);
            if (intersectioPoint.distanceFromCamera != INFINITY) {
                uint32_t colour = colourPalette(intersectioPoint.intersectedTriangle.colour);
                window.setPixelColour(x, y, colour);
            }
        }
    }
}

void rasterisedSceneWithShadow(std::vector<ModelTriangle> &modelTriangle,  std::vector<std::vector<float>> &depthBuffer, DrawingWindow &window) {
    for (int y = 0; y < HEIGHT; y++) {
        for (int x = 0; x < WIDTH; x++) {
            // converting from 2D pixel into 3D direction
            glm::vec3 rayDirection = pixelToDirection(x, y);
            RayTriangleIntersection closestIntersection = getClosestValidIntersection(modelTriangle, cameraPosition,rayDirection);

            glm::vec3 lightDirection = glm::normalize(closestIntersection.intersectionPoint - lightPosition);
            RayTriangleIntersection lightIntersection = getClosestValidIntersection(modelTriangle, lightPosition, lightDirection);

            if (closestIntersection.distanceFromCamera != INFINITY) {
                if (closestIntersection.triangleIndex == lightIntersection.triangleIndex) {
                    uint32_t colour = colourPalette(closestIntersection.intersectedTriangle.colour);
                    window.setPixelColour(x, y, colour);
                }
            }
        }
    }
}

void rasterisedSceneWithLighting(std::vector<ModelTriangle> &modelTriangle, std::vector<std::vector<float>> &depthBuffer, DrawingWindow &window) {
    for (int y = 0; y < HEIGHT; y++) {
        for (int x = 0; x < WIDTH; x++) {
            glm::vec3 rayDirection = pixelToDirection(x, y);
            RayTriangleIntersection closestIntersection = getClosestValidIntersection(modelTriangle, cameraPosition,rayDirection);
            Colour colour = closestIntersection.intersectedTriangle.colour;
            glm::vec3 normal = modelTriangleNormal(closestIntersection.intersectedTriangle);

            if (closestIntersection.distanceFromCamera != INFINITY) {
                glm::vec3 intersectionPoint = closestIntersection.intersectionPoint;
                glm::vec3 lightRay = lightPosition - intersectionPoint;
                float lightDistance = glm::length(lightRay);
                glm::vec3 lightDirection = glm::normalize(lightRay);
                glm::vec3 shadowRay = intersectionPoint + lightDirection * 0.001f;
                RayTriangleIntersection shadowIntersection = getClosestValidIntersection(modelTriangle, shadowRay, lightDirection);
                float shadowFactor = 1.0f;

                if (shadowIntersection.distanceFromCamera < lightDistance && shadowIntersection.triangleIndex != closestIntersection.triangleIndex) {
                    shadowFactor = 0.5f;
                }

                float proximityIntensity = proximityLighting(intersectionPoint);
                float diffuseIntensity = angleOfIncidenceLighting(intersectionPoint, normal);
                float specularIntensity = specularLighting(intersectionPoint, normal);

                glm::vec3 colourVec3 = glm::vec3(colour.red, colour.green, colour.blue);
                glm::vec3 lightingVec3 = (colourVec3 * (proximityIntensity * diffuseIntensity) + glm::vec3(255.0f) * specularIntensity) * shadowFactor;
                glm::vec3 ambient = colourVec3 * glm::vec3(ambientLighting);
                glm::vec3 lightingColour = glm::clamp(lightingVec3 + ambient, 0.0f, 255.0f);

                Colour combinedColour = Colour(lightingColour.x, lightingColour.y, lightingColour.z);
                window.setPixelColour(x, y, colourPalette(combinedColour));
            }
        }
    }
}

// Week7 Sphere

std::vector<ModelTriangle> readSphereOBJ(const std::string &fileName, const std::map<std::string, Colour> &readPalette, float scalingFactor) {
    std::vector<ModelTriangle> triangles;
    std::vector<glm::vec3> vertices;
    std::unordered_map<std::string, Colour> colours;
    std::ifstream file(fileName);
    std::string line;
    Colour colour = Colour(255, 0, 0);

    while (std::getline(file, line)) {
        std::vector<std::string> tokens = split(line, ' ');
        if (tokens.empty()) continue;
        if (tokens[0] == "v") {
            glm::vec3 vertex = glm::vec3(-scalingFactor * std::stof(tokens[1]),
                                         scalingFactor * std::stof(tokens[2]),
                                         scalingFactor * std::stof(tokens[3]));
            vertices.push_back(vertex);
        }
    }
    file.clear();
    file.seekg(0);
    while (std::getline(file, line)) {
        std::vector<std::string> tokens = split(line, ' ');
        if (tokens[0] == "f") {
            int v0 = std::stoi(split(tokens[1], '/')[0]) - 1;
            int v1 = std::stoi(split(tokens[2], '/')[0]) - 1;
            int v2 = std::stoi(split(tokens[3], '/')[0]) - 1;
            ModelTriangle triangle(vertices[v0], vertices[v1], vertices[v2], colour);
            triangle.normal = modelTriangleNormal(triangle);
            triangles.push_back(triangle);
        }
    }
    file.close();
    return triangles;
}

glm::vec3 vertexNormals(glm::vec3 &vertex, std::vector<ModelTriangle> &modelSphere) {
    glm::vec3 vertexNormal(0.0f, 0.0f, 0.0f);
    int count = 0;
    for (ModelTriangle &triangle : modelSphere) {
        if (std::find(triangle.vertices.begin(), triangle.vertices.end(), vertex) != triangle.vertices.end()) {
                vertexNormal += triangle.normal;
                count++;
        }
    }
    return count > 0 ? glm::normalize(vertexNormal / float(count)) : glm::vec3(0.0f, 0.0f, 0.0f);
}

Colour colourForLighting(Colour &colour, float lightIntensity, glm::vec3 &ambientLighting) {
    int red = std::min(roundToInt(colour.red * lightIntensity + 255 * ambientLighting.x), 255);
    int green = std::min(roundToInt(colour.green * lightIntensity + 255 * ambientLighting.y), 255);
    int blue = std::min(roundToInt(colour.blue * lightIntensity + 255 * ambientLighting.z), 255);
    return Colour(red, green, blue);
}

void sphereRasterisedSceneWithGourandShading(std::vector<ModelTriangle> &modelSphere, std::vector<std::vector<float>> &depthBuffer, DrawingWindow &window) {
    for (int y = 0; y < HEIGHT; y++) {
        for (int x = 0; x < WIDTH; x++) {
            glm::vec3 rayDirection = pixelToDirection(x, y);
            RayTriangleIntersection closestIntersection = getClosestValidIntersection(modelSphere, cameraPosition, rayDirection);

            if (closestIntersection.distanceFromCamera != INFINITY) {
                auto &triangle = closestIntersection.intersectedTriangle;
                Colour colour = triangle.colour;
                glm::vec3 intersectionPoint = closestIntersection.intersectionPoint;

                glm::vec3 v0 = triangle.vertices[0];
                glm::vec3 v1 = triangle.vertices[1];
                glm::vec3 v2 = triangle.vertices[2];
                glm::vec3 n0 = vertexNormals(v0, modelSphere);
                glm::vec3 n1 = vertexNormals(v1, modelSphere);
                glm::vec3 n2 = vertexNormals(v2, modelSphere);

                // Barycentric Coordinates
                float det = ((v1.y - v2.y) * (v0.x - v2.x)) + ((v2.x - v1.x) * (v0.y - v2.y));
                float u = (((v1.y - v2.y) * (intersectionPoint.x - v2.x)) + ((v2.x - v1.x) * (intersectionPoint.y - v2.y))) / det;
                float v = (((v2.y - v0.y) * (intersectionPoint.x - v2.x)) + ((v0.x - v2.x) * (intersectionPoint.y - v2.y))) / det;
                float w = 1 - u - v;
                glm::vec3 pointNormal = u * n0 + v * n1 + w * n2;

                float proximityIntensity = proximityLighting(intersectionPoint);
                float diffuseIntensity = angleOfIncidenceLighting(intersectionPoint, pointNormal);
                float specularIntensity = specularLighting(intersectionPoint, pointNormal);

                glm::vec3 colourVec3 = glm::vec3(colour.red, colour.green, colour.blue);
                glm::vec3 lightingVec3 = (colourVec3 * (proximityIntensity * diffuseIntensity) + glm::vec3(255.0f) * specularIntensity);
                glm::vec3 ambient = colourVec3 * glm::vec3(ambientLighting);
                glm::vec3 lightingColour = glm::clamp(lightingVec3 + ambient, 0.0f, 255.0f);

                Colour combinedColour = Colour(roundToInt(lightingColour.x), roundToInt(lightingColour.y), roundToInt(lightingColour.z));
                window.setPixelColour(x, y, colourPalette(combinedColour));
            }
        }
    }
}

void sphereRasterisedSceneWithFlat(std::vector<ModelTriangle> &modelSphere, std::vector<std::vector<float>> &depthBuffer, DrawingWindow &window) {
    for (int y = 0; y < HEIGHT; y++) {
        for (int x = 0; x < WIDTH; x++) {
            glm::vec3 rayDirection = pixelToDirection(x, y);
            RayTriangleIntersection closestIntersection = getClosestValidIntersection(modelSphere, cameraPosition,rayDirection);
            Colour colour = closestIntersection.intersectedTriangle.colour;
            glm::vec3 normal = modelTriangleNormal(closestIntersection.intersectedTriangle);

            if (closestIntersection.distanceFromCamera != INFINITY) {
                glm::vec3 intersectionPoint = closestIntersection.intersectionPoint;
                glm::vec3 lightRay = lightPosition - intersectionPoint;
                float lightDistance = glm::length(lightRay);
                glm::vec3 lightDirection = glm::normalize(lightRay);
                glm::vec3 shadowRay = intersectionPoint + lightDirection * 0.001f;
                RayTriangleIntersection shadowIntersection = getClosestValidIntersection(modelSphere, shadowRay, lightDirection);
                float shadowFactor = 1.0f;

                if (shadowIntersection.distanceFromCamera < lightDistance && shadowIntersection.triangleIndex != closestIntersection.triangleIndex) {
                    shadowFactor = 0.5f;
                }
                float proximityIntensity = proximityLighting(intersectionPoint);
                float diffuseIntensity = angleOfIncidenceLighting(intersectionPoint, normal);
                float specularIntensity = specularLighting(intersectionPoint, normal);

                glm::vec3 colourVec3 = glm::vec3(colour.red, colour.green, colour.blue);
                glm::vec3 lightingVec3 = (colourVec3 * (proximityIntensity * diffuseIntensity) + glm::vec3(255.0f) * specularIntensity) * shadowFactor;
                glm::vec3 ambient = colourVec3 * glm::vec3(ambientLighting);
                glm::vec3 lightingColour = glm::clamp(lightingVec3 + ambient, 0.0f, 255.0f);

                Colour combinedColour = Colour(lightingColour.x, lightingColour.y, lightingColour.z);
                window.setPixelColour(x, y, colourPalette(combinedColour));
            }
        }
    }
}

// Draw Controller

void drawWireframe(DrawingWindow &window) {
    window.clearPixels();
    std::vector<std::vector<float>> depthBuffer = initialiseDepthBuffer(WIDTH, HEIGHT);
    std::map<std::string, Colour> readPalette = readMTL("./src/cornell-box.mtl");
    std::vector<ModelTriangle> modelTriangle = readOBJ("./src/cornell-box.obj", readPalette, 0.35f);
    wireframeRender(modelTriangle, depthBuffer, window);
}

void drawRasterised(DrawingWindow &window) {
    window.clearPixels();
    std::vector<std::vector<float>> depthBuffer = initialiseDepthBuffer(WIDTH, HEIGHT);
    std::map<std::string, Colour> readPalette = readMTL("./src/cornell-box.mtl");
    std::vector<ModelTriangle> modelTriangle = readOBJ("./src/cornell-box.obj", readPalette, 0.35f);
    rasterisedRender(modelTriangle, depthBuffer, window);
    if (orbitAnimation) {
        orbitCamera();
    }
}

void drawRayScene(DrawingWindow &window) {
    window.clearPixels();
    std::vector<std::vector<float>> depthBuffer = initialiseDepthBuffer(WIDTH, HEIGHT);
    std::map<std::string, Colour> readPalette = readMTL("./src/cornell-box.mtl");
    std::vector<ModelTriangle> modelTriangle = readOBJ("./src/cornell-box.obj", readPalette, 0.35f);
    rasterisedScene(modelTriangle, depthBuffer, window);
}

void drawShadow(DrawingWindow &window) {
    window.clearPixels();
    std::vector<std::vector<float>> depthBuffer = initialiseDepthBuffer(WIDTH, HEIGHT);
    std::map<std::string, Colour> readPalette = readMTL("./src/cornell-box.mtl");
    std::vector<ModelTriangle> modelTriangle = readOBJ("./src/cornell-box.obj", readPalette, 0.35f);
    rasterisedSceneWithShadow(modelTriangle, depthBuffer, window);
}

void drawLighting(DrawingWindow &window) {
    window.clearPixels();
    std::vector<std::vector<float>> depthBuffer = initialiseDepthBuffer(WIDTH, HEIGHT);
    std::map<std::string, Colour> readPalette = readMTL("./src/cornell-box.mtl");
    std::vector<ModelTriangle> modelTriangle = readOBJ("./src/cornell-box.obj", readPalette, 0.35f);
    rasterisedSceneWithLighting(modelTriangle, depthBuffer, window);
}

void drawSphereWithFlat(DrawingWindow &window) {
    window.clearPixels();
    std::vector<std::vector<float>> depthBuffer = initialiseDepthBuffer(WIDTH, HEIGHT);
    std::map<std::string, Colour> readPalette = readMTL("./src/cornell-box.mtl");
    std::vector<ModelTriangle> readSphere = readSphereOBJ("./src/sphere.obj", readPalette, 0.50f);
    sphereRasterisedSceneWithFlat(readSphere, depthBuffer, window);
}

void drawSphereWithGourandShading (DrawingWindow &window) {
    window.clearPixels();
    std::vector<std::vector<float>> depthBuffer = initialiseDepthBuffer(WIDTH, HEIGHT);
    std::map<std::string, Colour> readPalette = readMTL("./src/cornell-box.mtl");
    std::vector<ModelTriangle> readSphere = readSphereOBJ("./src/sphere.obj", readPalette, 0.50f);
    sphereRasterisedSceneWithGourandShading(readSphere, depthBuffer, window);
}


void handleEvent(SDL_Event event, DrawingWindow &window) {
    if (event.type == SDL_KEYDOWN) {
        float lightMoveSpeed = 0.1f;
        float translationSpeed = 0.1f;
        float rotationSpeed = 2.0f;
        if (event.key.keysym.sym == SDLK_LEFT) {
            cameraTranslation(translationSpeed, 0.0f, 0.0f);
            std::cout << "LEFT" << std::endl;
        }
        else if (event.key.keysym.sym == SDLK_RIGHT) {
            cameraTranslation(-translationSpeed, 0.0f, 0.0f);
            std::cout << "RIGHT" << std::endl;
        }
        else if (event.key.keysym.sym == SDLK_UP) {
            cameraTranslation( 0.0f, translationSpeed, 0.0f);
            std::cout << "UP" << std::endl;
        }
        else if (event.key.keysym.sym == SDLK_DOWN) {
            cameraTranslation( 0.0f, -translationSpeed, 0.0f);
            std::cout << "DOWN" << std::endl;
        }
        else if (event.key.keysym.sym == SDLK_f) {
            cameraTranslation( 0.0f, 0.0f, -translationSpeed);
            std::cout << "FORWARDS" << std::endl;
        }
        else if (event.key.keysym.sym == SDLK_b) {
            cameraTranslation( 0.0f, 0.0f, translationSpeed);
            std::cout << "BACKWARDS" << std::endl;
        }
        else if (event.key.keysym.sym == SDLK_i) {
            cameraRotation(rotationSpeed, 0.0f);
        }
        else if (event.key.keysym.sym == SDLK_k) {
            cameraRotation(-rotationSpeed, 0.0f);
        }
        else if (event.key.keysym.sym == SDLK_j) {
            cameraRotation(0.0f, -rotationSpeed);
        }
        else if (event.key.keysym.sym == SDLK_l) {
            cameraRotation(0.0f, rotationSpeed);
        }
        else if (event.key.keysym.sym == SDLK_o) {
            orbitAnimation = !orbitAnimation;
            std::cout << "Orbit Camera " << (orbitAnimation ? "ON" : "OFF") << std::endl;
        }
        else if (event.key.keysym.sym == SDLK_a) {
            lightPosition.x -= lightMoveSpeed;
        }
        else if (event.key.keysym.sym == SDLK_d) {
            lightPosition.x += lightMoveSpeed;
        }
        else if (event.key.keysym.sym == SDLK_w) {
            lightPosition.y += lightMoveSpeed;
        }
        else if (event.key.keysym.sym == SDLK_s) {
            lightPosition.y -= lightMoveSpeed;
        }
        else if (event.key.keysym.sym == SDLK_q) {
            lightPosition.z -= lightMoveSpeed;
        }
        else if (event.key.keysym.sym == SDLK_e) {
            lightPosition.z += lightMoveSpeed;
        }
        else if (event.key.keysym.sym == SDLK_1) {
            callDraw = 1;
        }
        else if (event.key.keysym.sym == SDLK_2) {
            callDraw = 2;
        }
        else if (event.key.keysym.sym == SDLK_3) {
            callDraw = 3;
        }
        else if (event.key.keysym.sym == SDLK_4) {
            callDraw = 4;
        }
        else if (event.key.keysym.sym == SDLK_5) {
            callDraw = 5;
        }
    } else if (event.type == SDL_MOUSEBUTTONDOWN) {
        window.savePPM("output.ppm");
        window.saveBMP("output.bmp");
    }
}

int main(int argc, char *argv[]) {
    DrawingWindow window = DrawingWindow(WIDTH, HEIGHT, false);

    SDL_Event event;

    while (true) {
        // We MUST poll for events - otherwise the window will freeze !
        if (window.pollForInputEvents(event)) handleEvent(event, window);
        // drawSphereWithFlat(window);
        drawSphereWithGourandShading(window);

        /*
        if (callDraw == 1) {
            drawWireframe(window);
        }
        if (callDraw == 2) {
            drawRasterised(window);
        }
        if (callDraw == 3) {
             drawRayScene(window);
        }
        if (callDraw == 4) {
             drawShadow(window);
        }
        if (callDraw == 5) {
             drawLighting(window);
        }
        */

        // Need to render the frame at the end, or nothing actually gets shown on the screen !
        window.renderFrame();
    }
}
