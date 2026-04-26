#ifndef VECTOR_H
#define VECTOR_H

struct Vector2D {
    float x, y;

    Vector2D(float x = 0, float y = 0) : x(x), y(y) {}

    Vector2D operator+(const Vector2D& other) const {
        return {x + other.x, y + other.y};
    }

    Vector2D operator-(const Vector2D& other) const {
        return {x - other.x, y - other.y};
    }

    Vector2D operator*(float scalar) const {
        return {x * scalar, y * scalar};
    }

    Vector2D operator/(float scalar) const {
        return {x / scalar, y / scalar};
    }

    Vector2D operator+(float scalar) const {
        return {x + scalar, y + scalar};
    }

    Vector2D operator-(float scalar) const {
        return {x - scalar, y - scalar};
    }

    friend Vector2D operator*(float scalar, const Vector2D& vector) {
        return vector * scalar;
    }

    friend Vector2D operator+(float scalar, const Vector2D& vector) {
        return vector + scalar;
    }

    friend Vector2D operator-(float scalar, const Vector2D& vector) {
        return {scalar - vector.x, scalar - vector.y};
    }

    Vector2D& operator+=(const Vector2D& other) {
        x += other.x;
        y += other.y;
        return *this;
    }

    Vector2D& operator-=(const Vector2D& other) {
        x -= other.x;
        y -= other.y;
        return *this;
    }

    Vector2D& operator*=(float scalar) {
        x *= scalar;
        y *= scalar;
        return *this;
    }

    Vector2D& operator/=(float scalar) {
        x /= scalar;
        y /= scalar;
        return *this;
    }

    Vector2D& operator+=(float scalar) {
        x += scalar;
        y += scalar;
        return *this;
    }

    Vector2D& operator-=(float scalar) {
        x -= scalar;
        y -= scalar;
        return *this;
    }
};

#endif
