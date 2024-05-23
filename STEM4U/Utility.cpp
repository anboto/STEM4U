// SPDX-License-Identifier: Apache-2.0
// Copyright 2021 - 2022, the Anboto author and contributors
#include <Core/Core.h>
#include <Functions4U/Functions4U.h>
#include <Eigen/Eigen.h>
#include "Utility.h"

namespace Upp {

using namespace Eigen;


Image ApplyHomography(const Image& orig, const Color &back, 
			const Point &from0, const Point &from1, const Point &from2, const Point &from3,
			const Point &to0,   const Point &to1,   const Point &to2,   const Point &to3) {
	Homography<int, double> h(to0, to1, to2, to3, from0, from1, from2, from3);
	
	Rect rcdest(min(to0.x, to1.x, to2.x, to3.x), min(to0.y, to1.y, to2.y, to3.y), max(to0.x, to1.x, to2.x, to3.x), max(to0.y, to1.y, to2.y, to3.y));
	Size szdest(rcdest.Width(), rcdest.Height());
	
	ImageBuffer ib(szdest);
	ib.SetKind(IMAGE_OPAQUE);
	if (!IsNull(back)) {
		for(int y = 0; y < szdest.cy; y++)
			Fill(ib[y], back, szdest.cx);	
	}
	
	int dety0 = 0;
	for(int desty = rcdest.top; desty < rcdest.bottom; desty++, dety0++) {
		int destx0 = 0;
		for(int destx = rcdest.left; destx < rcdest.right; destx++, destx0++) {
			Point from = h.Transform(Point(destx, desty));
			if (from.x >= 0 && from.y >= 0 && from.x < orig.GetWidth() && from.y < orig.GetHeight()) {
				RGBA *s = ib[dety0] + destx0;
				*s = *GetPixel(orig, from.x, from.y);
			}
		}
	}
	return ib;
}

}

