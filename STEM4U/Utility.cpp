// SPDX-License-Identifier: Apache-2.0
// Copyright 2021 - 2022, the Anboto author and contributors
#include <Core/Core.h>
#include <Functions4U/Functions4U.h>
#include <Eigen/Eigen.h>
#include <ScatterDraw/DataSource.h>
#include "Utility.h"

#include <Functions4U/EnableWarnings.h>

namespace Upp {
using namespace Eigen;

RGBA GetPixelBilinear(const Image &img, double x, double y) {
	int x1 = int(x);
	int x2;
	if (x1+1 == img.GetWidth())
		x2 = x1;
	else
		x2 = x1 + 1;
	int y1 = int(y);
	int y2;
	if (y1+1 == img.GetHeight())
		y2 = y1;
	else
		y2 = y1 + 1;

	byte z11, z12, z21, z22;
	RGBA ret;
	
	z11 = img[y1][x1].r;
	z12 = img[y2][x1].r;
	z21 = img[y1][x2].r;
	z22 = img[y2][x2].r;
	ret.r = BilinearInterpolate<double, int>(x, y, x1, x2, y1, y2, z11, z12, z21, z22);
	
	z11 = img[y1][x1].g;
	z12 = img[y2][x1].g;
	z21 = img[y1][x2].g;
	z22 = img[y2][x2].g;
	ret.g = BilinearInterpolate<double, int>(x, y, x1, x2, y1, y2, z11, z12, z21, z22);
	
	z11 = img[y1][x1].b;
	z12 = img[y2][x1].b;
	z21 = img[y1][x2].b;
	z22 = img[y2][x2].b;
	ret.b = BilinearInterpolate<double, int>(x, y, x1, x2, y1, y2, z11, z12, z21, z22);
	
	ret.a = 255;
	
	return ret;
}

RGBA GetPixel(const Image &img, double x, double y, bool bilinear) {
	if (bilinear)
		return GetPixelBilinear(img, x, y);
	else
		return RGBA(GetPixelBilinear(img, (int)x, (int)y));
}
	
Image ApplyHomography(const Image& orig, const Color &back, 
			const Point &from0, const Point &from1, const Point &from2, const Point &from3,
			const Point &to0,   const Point &to1,   const Point &to2,   const Point &to3, bool bilinear, bool fit) {
	Homography<int> h;
	
	h.QuadToQuad(to0, to1, to2, to3, from0, from1, from2, from3);
	
	Rect rcdest;
	if (fit)
		rcdest = Rect(min(to0.x, to1.x, to2.x, to3.x), min(to0.y, to1.y, to2.y, to3.y), max(to0.x, to1.x, to2.x, to3.x), max(to0.y, to1.y, to2.y, to3.y));
	else
		rcdest = Rect(0, 0, max(to0.x, to1.x, to2.x, to3.x), max(to0.y, to1.y, to2.y, to3.y));
	
	Size szdest(rcdest.Width(), rcdest.Height());
	
	ImageBuffer ib(szdest);
	ib.SetKind(IMAGE_OPAQUE);
	if (!IsNull(back)) {
		for(int y = 0; y < szdest.cy; y++)
			Fill(ib[y], back, (size_t)szdest.cx);	
	}
	
	int dety0 = 0;
	for(int desty = rcdest.top; desty < rcdest.bottom; desty++, dety0++) {
		int destx0 = 0;
		for(int destx = rcdest.left; destx < rcdest.right; destx++, destx0++) {
			Pointf from = h.Transform(Point(destx, desty));
			if (from.x >= 0 && from.y >= 0 && from.x < orig.GetWidth() && from.y < orig.GetHeight()) {
				RGBA *s = ib[dety0] + destx0;
				*s = GetPixel(orig, from.x, from.y, bilinear);
			}
		}
	}
	return ib;
}

Image ApplyHomography(const Image& orig, const Color &back, 
			const Point &from0, const Point &from1, const Point &from2, const Point &from3, const Size &sz, bool bilinear) {
	return ApplyHomography(orig, back, from0, from1, from2, from3, Point(0, 0), Point(sz.cx, 0), Point(sz.cx, sz.cy), Point(0, sz.cy), bilinear);
}

static void MapToSphere(const Point &p, const Size &sz, Vector3d& v) {	// Maps 2D mouse to a virtual unit sphere (trackball), centered in the screen.
    double nx = (sz.cx - 2 * p.x)/(double)sz.cx;// Normalise x and y to range [-1, 1], with origin at screen centre
    double ny = (sz.cy - 2 * p.y)/(double)sz.cy;

    double len2 = nx*nx + ny*ny;				// Compute squared distance from centre

    if (len2 > 1) {
        double norm = 1/::sqrt(len2);			// Outside unit circle: project onto the sphere by normalisation
        v = Vector3d(nx*norm, ny*norm, 0);		// lies on circle in XY plane
    } else
        v = Vector3d(nx, ny, ::sqrt(1 - len2));	// Inside unit circle: point on hemisphere (z > 0)
}

Affine3d TrackballRotation(const Point &p0, const Point &p1, const Size &sz, const Vector3d& centre) {
    Vector3d v0, v1;
    MapToSphere(p0, sz, v0);
    MapToSphere(p1, sz, v1);

    double angle = ::acos(std::clamp(v0.dot(v1), -1., 1.));	// Get the angle between the two vectors using the dot product
	
	Vector3d axis = v0.cross(v1);							// Get the rotation axis using the cross product
	axis.normalize();

	Affine3d ret = Affine3d::Identity();
	
    if (axis.norm() < 1e-6 || ::abs(angle) < 1e-6)			// No rotation
        return ret;

    return ret.translate(centre).rotate(AngleAxisd(angle, axis)).translate(-centre);
}



// Detect grid structures with tolerance in spacing
bool DetectGrid(UVector<Pointf>& pts, double tol, Pointf &topLeft, Pointf &bottomRight, int &cols, int &rows, UVector<int> &ids) {
	auto FindWithTolerance = [](const UVector<Pointf>& idx, const Pointf& p, double tol) {
	    for (int i = 0; i < idx.size(); i++) {
	        if (abs(idx[i].x - p.x) <= tol && abs(idx[i].y - p.y) <= tol)
	            return i;
	    }
	    return -1;
	};
	struct GridFit {
	    Pointf topLeft, bottomRight;
	    int cols, rows;
	    UVector<int> ids; // The points forming this grid
	};
	
    GridFit best;
    int bestCount = 0;

    // Sort input for easier grouping
    Sort(pts, [](const Pointf& a, const Pointf& b) {
        if(a.y != b.y) 
        	return a.y < b.y;
        return a.x < b.x;
    });

	bool found = false;
	
    // Try each point as potential top-left corner
    for (int i = 0; i < pts.size() && !found; i++) {
        for (int j = i+1; j < pts.size() && !found; j++) {
            if (abs(pts[i].y - pts[j].y) > tol) 
            	break; // must be same row within tolerance
            
            double dx = pts[j].x - pts[i].x;
            if (dx <= 0) 
            	continue;

            // Try to detect vertical spacing
            for (int k = i+1; k < pts.size(); k++) {
                if (abs(pts[k].x - pts[i].x) <= tol && pts[k].y > pts[i].y) {
                    double dy = pts[k].y - pts[i].y;
                    if (dy <= 0) 
                    	continue;

                    GridFit g;				// Try to expand a grid with spacing dx, dy
                    g.topLeft = pts[i];
                    g.cols = 0;
                    g.rows = 0;

                    UVector<Pointf> idx = clone(pts);

                    int col = 0;
                    while(FindWithTolerance(idx, Pointf(g.topLeft.x + col*dx, g.topLeft.y), tol) >= 0) 
                    	col++;
                    int row = 0;
                    while(FindWithTolerance(idx, Pointf(g.topLeft.x, g.topLeft.y + row*dy), tol) >= 0) 
                    	row++;

                    g.cols = col;
                    g.rows = row;
					
                    for (int r = 0; r < row; r++) {			// Collect valid grid points
                        for (int c = 0; c < col; c++) {
                            Pointf test(g.topLeft.x + c*dx, g.topLeft.y + r*dy);
                            int id;
                            if((id = FindWithTolerance(idx, test, tol)) >= 0)
                                g.ids << id;
                            else {		// Break, there is a hole in the grid
                                r = row;
                                c = col;
                                g.ids.Clear();
                            }
                        }
                    }
                    if(g.ids.size() > bestCount) {
                        g.bottomRight = pick(idx[FindWithTolerance(idx, Pointf(g.topLeft.x + (col-1)*dx, g.topLeft.y + (row-1)*dy), tol)]);
                        best = pick(g);
                        bestCount = best.ids.size();
                        if (pts.size() <= best.rows*(best.cols+1) || pts.size() <= (best.rows+1)*best.cols) {	// Impossible to get a bigger grid
                            found = true;
                            break;
                        }
                    }
                }
            }
        }
    }
    if (bestCount == 0)
        return false;
    
    Sort(best.ids);
    
    topLeft = best.topLeft;
    bottomRight = best.bottomRight;
    cols = best.cols;
    rows = best.rows;
    ids = pick(best.ids);
    
    return true;
}

}
