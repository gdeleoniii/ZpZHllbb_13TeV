int getIntersection(TGraph *thisGraph, TGraph *thatGraph, vector<double> *intersectX, vector<double> *intersectY){

  int intersectionPt = 0;
  
  // Loop over all points in thisGraph
  for( int gi = 0; gi < thisGraph->GetN()-1; ++gi ){
    
    // Loop over all points in thatGraph
    for( int gj = 0; gj < thatGraph->GetN()-1; ++gj ){

      // Get the current point and the next point for each graph
      double ix1, iy1, ix2, iy2;
      double jx1, jy1, jx2, jy2;

      thisGraph->GetPoint(gi,   ix1, iy1);
      thisGraph->GetPoint(gi+1, ix2, iy2);
      thatGraph->GetPoint(gj,   jx1, jy1);
      thatGraph->GetPoint(gj+1, jx2, jy2);

      // Calculate the intersection between two straight lines, x axis
      double myX = (jx1*(jy2*(ix1-ix2)+ix2*iy1-ix1*iy2)+jx2*(jy1*(-ix1+ix2)-ix2*iy1+ix1*iy2))/(-(jy1-jy2)*(ix1-ix2)+(jx1-jx2)*(iy1-iy2));

      // Calculate the intersection between two straight lines, y axis
      double myY = (jx1*jy2*(iy1-iy2)+jx2*jy1*(-iy1+iy2)+(jy1-jy2)*(ix2*iy1-ix1*iy2))/(-(jy1-jy2)*(ix1-ix2)+(jx1-jx2)*(iy1-iy2));

      // Find the tightest interval along the x-axis defined by the four points
      double xMin = max(min(ix1, ix2), min(jx1, jx2));
      double xMax = min(max(ix1, ix2), max(jx1, jx2));

      // If points from the two lines overlap, they are trivially intersecting
      if( (ix1 == jx1 && iy1 == jy1) || (ix2 == jx2 && iy2 == jy2) ){

	intersectX->push_back((ix1==jx1 && iy1==jy1) ? ix1 : ix2);
	intersectY->push_back((ix1==jx1 && iy1==jy1) ? iy1 : iy2);

	++intersectionPt;
	
      }

      // If the intersection between the two lines is within the tight range, add it to the list of intersections.
      else if( myX > xMin && myX < xMax ){

	intersectX->push_back(myX);
	intersectY->push_back(myY);

	++intersectionPt;
	
      }
      
    } // end of thatGraph
    
  } // end of thisGraph

  // Return number of intersection points
  return intersectionPt;
  
}
