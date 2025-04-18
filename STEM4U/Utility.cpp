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

}
