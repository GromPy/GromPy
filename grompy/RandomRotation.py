#Taken from graphicgems (http://www.graphicsgems.org/)
#
#
#GraphiGems.h:
#/*
# * GraphicsGems.h
# * Version 1.0 - Andrew Glassner
# * from "Graphics Gems", Academic Press, 1990
# */
#
##ifndef GG_H
#
##define GG_H 1
#
#/*********************/
#/* 2d geometry types */
#/*********************/
#
#typedef struct Point2Struct {    /* 2d point */
#    double x, y;
#    } Point2;
#typedef Point2 Vector2;
#
#typedef struct IntPoint2Struct {    /* 2d integer point */
#    int x, y;
#    } IntPoint2;
#
#typedef struct Matrix3Struct {    /* 3-by-3 matrix */
#    double element[3][3];
#    } Matrix3;
#
#typedef struct Box2dStruct {        /* 2d box */
#    Point2 min, max;
#    } Box2;
#
#
#/*********************/
#/* 3d geometry types */
#/*********************/
#
#typedef struct Point3Struct {    /* 3d point */
#    double x, y, z;
#    } Point3;
#typedef Point3 Vector3;
#
#typedef struct IntPoint3Struct {    /* 3d integer point */
#    int x, y, z;
#    } IntPoint3;
#
#
#typedef struct Matrix4Struct {    /* 4-by-4 matrix */
#    double element[4][4];
#    } Matrix4;
#
#typedef struct Box3dStruct {        /* 3d box */
#    Point3 min, max;
#    } Box3;
#
#
#
#/***********************/
#/* one-argument macros */
#/***********************/
#
#/* absolute value of a */
##define ABS(a)        (((a)<0) ? -(a) : (a))
#
#/* round a to nearest int */
##define ROUND(a)    floor((a)+0.5)
#
#/* take sign of a, either -1, 0, or 1 */
##define ZSGN(a)        (((a)<0) ? -1 : (a)>0 ? 1 : 0)
#
#/* take binary sign of a, either -1, or 1 if >= 0 */
##define SGN(a)        (((a)<0) ? -1 : 1)
#
#/* shout if something that should be true isn't */
##define ASSERT(x) \
#if (!(x)) fprintf(stderr," Assert failed: x\n");
#
#/* square a */
##define SQR(a)        ((a)*(a))
#
#
#/***********************/
#/* two-argument macros */
#/***********************/
#
#/* find minimum of a and b */
##define MIN(a,b)    (((a)<(b))?(a):(b))
#
#/* find maximum of a and b */
##define MAX(a,b)    (((a)>(b))?(a):(b))
#
#/* swap a and b (see Gem by Wyvill) */
##define SWAP(a,b)    { a^=b; b^=a; a^=b; }
#
#/* linear interpolation from l (when a=0) to h (when a=1)*/
#/* (equal to (a*h)+((1-a)*l) */
##define LERP(a,l,h)    ((l)+(((h)-(l))*(a)))
#
#/* clamp the input to the specified range */
##define CLAMP(v,l,h)    ((v)<(l) ? (l) : (v) > (h) ? (h) : v)
#
#
#/****************************/
#/* memory allocation macros */
#/****************************/
#
#/* create a new instance of a structure (see Gem by Hultquist) */
##define NEWSTRUCT(x)    (struct x *)(malloc((unsigned)sizeof(struct x)))
#
#/* create a new instance of a type */
##define NEWTYPE(x)    (x *)(malloc((unsigned)sizeof(x)))
#
#
#/********************/
#/* useful constants */
#/********************/
#
##define PI        3.141592    /* the venerable pi */
##define PITIMES2    6.283185    /* 2 * pi */
##define PIOVER2        1.570796    /* pi / 2 */
##define E        2.718282    /* the venerable e */
##define SQRT2        1.414214    /* sqrt(2) */
##define SQRT3        1.732051    /* sqrt(3) */
##define GOLDEN        1.618034    /* the golden ratio */
##define DTOR        0.017453    /* convert degrees to radians */
##define RTOD        57.29578    /* convert radians to degrees */
#
#
#/************/
#/* booleans */
#/************/
#
##define TRUE        1
##define FALSE        0
##define ON        1
##define OFF         0
#typedef int boolean;            /* boolean data type */
#typedef boolean flag;            /* flag data type */
#
#extern double V2SquaredLength(), V2Length();
#extern double V2Dot(), V2DistanceBetween2Points();
#extern Vector2 *V2Negate(), *V2Normalize(), *V2Scale(), *V2Add(), *V2Sub();
#extern Vector2 *V2Lerp(), *V2Combine(), *V2Mul(), *V2MakePerpendicular();
#extern Vector2 *V2New(), *V2Duplicate();
#extern Point2 *V2MulPointByProjMatrix();
#extern Matrix3 *V2MatMul(), *TransposeMatrix3();
#
#extern double V3SquaredLength(), V3Length();
#extern double V3Dot(), V3DistanceBetween2Points();
#extern Vector3 *V3Normalize(), *V3Scale(), *V3Add(), *V3Sub();
#extern Vector3 *V3Lerp(), *V3Combine(), *V3Mul(), *V3Cross();
#extern Vector3 *V3New(), *V3Duplicate();
#extern Point3 *V3MulPointByMatrix(), *V3MulPointByProjMatrix();
#extern Matrix4 *V3MatMul();
#
#extern double RegulaFalsi(), NewtonRaphson(), findroot();
#
##endif
#
#
#rotate.c:
##include <math.h>
##include "GraphicsGems.h"
#/*======================================================================*
# *  R A N D _ R O T A T I O N      Author: Jim Arvo, 1991               *
# *                                                                      *
# *  This routine maps three values (x[0], x[1], x[2]) in the range      *
# *  [0,1] into a 3x3 rotation matrix M.  Uniformly distributed random   *
# *  variables x0, x1, and x2 create uniformly distributed random        *
# *  rotation matrices.  To create small uniformly distributed           *
# *  "perturbations", supply samples in the following ranges             *
# *                                                                      *
# *      x[0] in [ 0, d ]                                                *
# *      x[1] in [ 0, 1 ]                                                *
# *      x[2] in [ 0, d ]                                                *
# *                                                                      *
# * where 0 < d < 1 controls the size of the perturbation.  Any of the   *
# * random variables may be stratified (or "jittered") for a slightly    *
# * more even distribution.                                              *
# *                                                                      *
# *======================================================================*/
#void rand_rotation( float x[], Matrix3 *M )
#    {
#    float theta = x[0] * PITIMES2; /* Rotation about the pole (Z).      */
#    float phi   = x[1] * PITIMES2; /* For direction of pole deflection. */
#    float z     = x[2] * 2.0;      /* For magnitude of pole deflection. */
#
#    /* Compute a vector V used for distributing points over the sphere  */
#    /* via the reflection I - V Transpose(V).  This formulation of V    */
#    /* will guarantee that if x[1] and x[2] are uniformly distributed,  */
#    /* the reflected points will be uniform on the sphere.  Note that V */
#    /* has length sqrt(2) to eliminate the 2 in the Householder matrix. */
#
#    float r  = sqrt( z );
#    float Vx = sin( phi ) * r;
#    float Vy = cos( phi ) * r;
#    float Vz = sqrt( 2.0 - z );
#
#    /* Compute the row vector S = Transpose(V) * R, where R is a simple */
#    /* rotation by theta about the z-axis.  No need to compute Sz since */
#    /* it's just Vz.                                                    */
#
#    float st = sin( theta );
#    float ct = cos( theta );
#    float Sx = Vx * ct - Vy * st;
#    float Sy = Vx * st + Vy * ct;
#
#    /* Construct the rotation matrix  ( V Transpose(V) - I ) R, which   */
#    /* is equivalent to V S - R.                                        */
#
#    M->element[0][0] = Vx * Sx - ct;
#    M->element[0][1] = Vx * Sy - st;
#    M->element[0][2] = Vx * Vz;
#
#    M->element[1][0] = Vy * Sx + st;
#    M->element[1][1] = Vy * Sy - ct;
#    M->element[1][2] = Vy * Vz;
#
#    M->element[2][0] = Vz * Sx;
#    M->element[2][1] = Vz * Sy;
#    M->element[2][2] = 1.0 - z;   /* This equals Vz * Vz - 1.0 */
#    }


import math
from scipy import array,rand,empty

PI       = 4.0 * math.atan(1.0)
PITIMES2 = 2.0*PI

def RandRotation(x = array):
    theta = x[0] * PITIMES2
    phi   = x[1] * PITIMES2
    z     = x[2] * 2.0

    r  = math.sqrt(z)
    Vx = math.sin(phi)*r
    Vy = math.cos(phi)*r
    Vz = math.sqrt(2.0-z)

    st = math.sin(theta)
    ct = math.cos(theta)
    Sx = Vx*ct - Vy*st
    Sy = Vx*st + Vy*ct

    M = empty((3,3))
    M[:] = 0.0

    M[0][0] = Vx*Sx - ct
    M[0][1] = Vx*Sy - st
    M[0][2] = Vx*Vz

    M[1][0] = Vy*Sx + st
    M[1][1] = Vy*Sy - ct
    M[1][2] = Vy*Vz

    M[2][0] = Vz*Sx
    M[2][1] = Vz*Sy
    M[2][2] = 1.0 - z

    return M

if __name__ == "__main__":
    "Do the work"

    print '*** hello world! ***'
#    a = array([random.random(),random.random(),random.random()])
    a = rand(3)
    print " a = "+str(a)+"\n"
    print "-a = "+str(-a)+"\n"
    print "RotMat = "+str(RandRotation(a))+"\n"


