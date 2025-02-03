#include <math.h>
#include <string>

struct Vector {
    float x;
    float y;
    float z;
};

struct Atom {
    const char* a;
    float x;
    float y;
    float z;
    float e;
    float s;
    float q;
};

Vector getPos(Atom a) {
    Vector output;
    output.x= a.x;
    output.y= a.y;
    output.z= a.z;
    return output;
}

Vector addVectors(Vector v1, Vector v2) {
    Vector output;
    output.x= v1.x + v2.x;
    output.y= v1.y + v2.y;
    output.z= v1.z + v2.z;
    return output;
}

Vector substracVectors(Vector v1, Vector v2) {
    Vector output;
    output.x= v1.x - v2.x;
    output.y= v1.y - v2.y;
    output.z= v1.z - v2.z;
    return output;
}

float magnitudeOfVector(Vector v) {
    return( sqrt(v.x*v.x + v.y*v.y + v.z*v.z) );
}

Vector productByScalar(Vector v, float k) {
    Vector output;
    output.x= v.x * k;
    output.y= v.y * k;
    output.z= v.z * k;
    return output;
}

Vector crossProduct(Vector v1, Vector v2) {
    Vector output;
    output.x= v1.y*v2.z-v1.z*v2.y;
    output.y= v1.z*v2.x-v1.x*v2.z;
    output.z= v1.x*v2.y-v1.y*v2.x;
    return output;
}

float distanceBetween(Vector c1, Vector c2, Vector bounds) {
    float dx= c1.x - c2.x;
    float dy= c1.y - c2.y;
    float dz= c1.z - c2.z;
    dx-= round(dx / bounds.x) * bounds.x;
    dy-= round(dy / bounds.y) * bounds.y;
    dz-= round(dz / bounds.z) * bounds.z;
    return sqrt(dx*dx + dy*dy + dz*dz);
}

float distanceBetween(Atom a, Vector v2, Vector bounds) {
    Vector v;
    v.x= a.x; v.y= a.y; v.z= a.z;
    return distanceBetween(v,v2,bounds);
}

float distanceBetween(Atom a1, Atom a2, Vector bounds) {
    Vector v;
    v.x= a2.x; v.y= a2.y; v.z= a2.z;
    return distanceBetween(a1,v,bounds);
}

float getAngle(Vector c1, Vector c2, Vector c3, Vector bounds) {
    float a= distanceBetween(c1,c3,bounds); //Opposite to the angle
    float b= distanceBetween(c1,c2,bounds);
    float c= distanceBetween(c2,c3,bounds);
    return abs(acos((pow(b,2)+pow(c,2)-pow(a,2))/(2*b*c))); //Cosine rule
}
