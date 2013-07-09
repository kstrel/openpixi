//Particle attributes indexes
#define P_SIZE 22

#define X 0
#define Y 1
#define Radius 2
#define Vx 3
#define Vy 4
#define Ax 5
#define Ay 6
#define Mass 7
#define Charge 8
#define PrevX 9
#define PrevY 10
#define Ex 11
#define Ey 12
#define PBz 13
#define PrevPositionComponentForceX 14
#define PrevPositionComponentForceY 15
#define PrevTangentVelocityComponentOfForceX 16
#define PrevTangentVelocityComponentOfForceY 17
#define PrevNormalVelocityComponentOfForceX 18
#define PrevNormalVelocityComponentOfForceY 19
#define PrevBz 20
#define PrevLinearDragCoefficient 21


//Force attributes indexes
#define ForceX 0
#define ForceY 1
#define PositionComponentofForceX 2
#define PositionComponentofForceY 3
#define TangentVelocityComponentOfForceX 4
#define TangentVelocityComponentOfForceY 5
#define NormalVelocityComponentofForceX 6
#define NormalVelocityComponentofForceY 7
#define Bz 8
#define LinearDragCoefficient 9

//boundary
#define X_MIN 0
#define X_CENTER 1
#define X_MAX 2
#define Y_MIN 0
#define Y_CENTER 3
#define Y_MAX 6

//Grid details
#define EXTRA_CELLS_BEFORE_GRID 2
#define EXTRA_CELLS_AFTER_GRID 2

//Cell attributes indexes
#define Cjx 0
#define Cjy 1
#define Crho 2
#define Cphi 3
#define Cex 4
#define Cey 5
#define Cbz 6
#define Cbzo 7


int get_region(double xmin, double xmax, double ymin, double ymax, double width, double height) {
    int xidx;
    int yidx;

    if (xmin < 0) {
            xidx  = X_MIN;
    } else if (xmax >= width) {
            xidx = X_MAX;
    } else {
            xidx = X_CENTER;
    }

    if (ymin < 0) {
            yidx = Y_MIN;
    } else if (ymax >= height) {
            yidx = Y_MAX;
    } else {
            yidx = Y_CENTER;
    }

    return xidx + yidx;
}

void fourBoundaryMove(int lx, int ly, double x, double y, double deltaX, 
                      double deltaY, double pCharge, double tstep, double *JxJy) {

    JxJy[0] = JxJy[0] + (lx, ly - 1, pCharge() * deltaX * ((1 - deltaY) / 2 - y));
    JxJy[0] = JxJy[0] + (lx, ly, p.getCharge() * deltaX * ((1 + deltaY) / 2 + y));
    JxJy[1] = JxJy[1] + (lx - 1, ly, p.getCharge() * deltaY * ((1 - deltaX) / 2 - x));
    JxJy[1] = JxJy[1] + (lx, ly, p.getCharge() * deltaY * ((1 + deltaX) / 2 + x));

}

private void sevenBoundaryMove(double x, double y, int xStart, int yStart, int xEnd, int yEnd,
                               double deltaX, double deltaY, double p, double tstep, double *JxJy) {
    //7-boundary move with equal y?
    if (yStart == yEnd) {
            //particle moves right?
            if (xEnd > xStart) {

                    double deltaX1 = 0.5 - x;
                    double deltaY1 = (deltaY / deltaX) * deltaX1;
                    fourBoundaryMove(xStart, yStart, x, y, deltaX1, deltaY1, p, tstep, JxJy);

                    deltaX -= deltaX1;
                    deltaY -= deltaY1;
                    y += deltaY1;
                    fourBoundaryMove(xEnd, yEnd, - 0.5, y, deltaX, deltaY, p, tstep, JxJy);

            }
            //particle moves left
            else {

                    double deltaX1 = -(0.5 + x);
                    double deltaY1 = (deltaY / deltaX) * deltaX1;
                    fourBoundaryMove(xStart, yStart, x, y, deltaX1, deltaY1, p, tstep, JxJy);

                    deltaX -= deltaX1;
                    deltaY -= deltaY1;
                    y += deltaY1;
                    fourBoundaryMove(xEnd, yEnd, 0.5, y, deltaX, deltaY, p, tstep, JxJy);

            }
    }
    //7-boundary move with equal x?
    if (xStart == xEnd) {
            //particle moves up?
            if (yEnd > yStart) {

                    double deltaY1 = 0.5 - y;
                    double deltaX1 = deltaX  * (deltaY1 / deltaY);
                    fourBoundaryMove(xStart, yStart, x, y, deltaX1, deltaY1, p, tstep, JxJy);

                    deltaX -= deltaX1;
                    deltaY -= deltaY1;
                    y += deltaY1;
                    fourBoundaryMove(xEnd, yEnd, x, -0.5, deltaX, deltaY, p, tstep, JxJy);

            }
            //particle moves down
            else {

                    double deltaY1 = -(0.5 + y);
                    double deltaX1 = (deltaX / deltaY) * deltaY1;
                    fourBoundaryMove(xStart, yStart, x, y, deltaX1, deltaY1, p, tstep, JxJy);

                    deltaX -= deltaX1;
                    deltaY -= deltaY1;
                    y += deltaY1;
                    fourBoundaryMove(xEnd, yEnd, x, 0.5, deltaX, deltaY, p, tstep, JxJy);

            }
    }

}

private void tenBoundaryMove(double x, double y, int xStart, int yStart, int xEnd, int yEnd,
			     double deltaX, double deltaY, double p, double tstep, double *JxJy) {
    //moved right?
    if (xEnd == (xStart+1)) {
            //moved up?
            if (yEnd == (yStart+1)) {

                    double deltaX1 = 0.5 - x;

                    //lower local origin
                    if(((deltaY / deltaX) * deltaX1 + y) < 0.5) {

                            double deltaY1 = (deltaY / deltaX) * deltaX1;
                            fourBoundaryMove(xStart, yStart, x, y, deltaX1, deltaY1, p, tstep, JxJy);

                            double deltaY2 = 0.5 - y - deltaY1;
                            double deltaX2 = (deltaX1 / deltaY1) * deltaY2;
                            y += deltaY1;
                            fourBoundaryMove(xStart+1, yStart, -0.5, y, deltaX2, deltaY2, p, tstep, JxJy);

                            deltaX -= (deltaX1 + deltaX2);
                            deltaY -= (deltaY1 + deltaY2);
                            x = deltaX2 - 0.5;
                            fourBoundaryMove(xEnd, yEnd, x, -0.5, deltaX, deltaY, p, tstep, JxJy);

                    }
                    //upper local origin
                    else {

                            double deltaY1 = 0.5 - y;
                            deltaX1 = (deltaX / deltaY) * deltaY1;
                            fourBoundaryMove(xStart, yStart, x, y, deltaX1, deltaY1, p, tstep, JxJy);

                            double deltaX2 = 0.5 - x - deltaX1;
                            double deltaY2 = (deltaY1 / deltaX1) * deltaX2;
                            x += deltaX1;
                            fourBoundaryMove(xStart, yStart+1, x, -0.5, deltaX2, deltaY2, p, tstep, JxJy);

                            deltaX -= (deltaX1 + deltaX2);
                            deltaY -= (deltaY1 + deltaY2);
                            y = deltaY2 - 0.5;
                            fourBoundaryMove(xEnd, yEnd, -0.5, y, deltaX, deltaY, p, tstep, JxJy);

                    }
            }
            //moved down
            else {

                    double deltaY1 = -(0.5 + y);

                    //lower local origin
                    if(((deltaX / deltaY) * deltaY1 + x) < 0.5) {

                            double deltaX1 = (deltaX / deltaY) * deltaY1;
                            fourBoundaryMove(xStart, yStart, x, y, deltaX1, deltaY1, p, tstep, JxJy);

                            double deltaX2 = 0.5 - x - deltaX1;
                            double deltaY2 = (deltaY / deltaX) * deltaX2;
                            x += deltaX1;
                            fourBoundaryMove(xStart, yStart-1, x, 0.5, deltaX2, deltaY2, p, tstep, JxJy);

                            deltaX -= (deltaX1 + deltaX2);
                            deltaY -= (deltaY1 + deltaY2);
                            y = 0.5 + deltaY2;
                            fourBoundaryMove(xEnd, yEnd, -0.5, y, deltaX, deltaY, p, tstep, JxJy);

                    }
                    //upper local origin
                    else {

                            double deltaX1 = 0.5 - x;
                            deltaY1 = (deltaY / deltaX) * deltaX1;
                            fourBoundaryMove(xStart, yStart, x, y, deltaX1, deltaY1, p, tstep, JxJy);

                            double deltaY2 = -(0.5 + y + deltaY1);
                            double deltaX2 = (deltaX1 / deltaY1) * deltaY2;
                            y += deltaY1;
                            fourBoundaryMove(xStart+1, yStart, -0.5, y, deltaX2, deltaY2, p, tstep, JxJy);

                            deltaX -= (deltaX1 + deltaX2);
                            deltaY -= (deltaY1 + deltaY2);
                            x = deltaX2 - 0.5;
                            fourBoundaryMove(xEnd, yEnd, x, 0.5, deltaX, deltaY, p, tstep, JxJy);

                    }
            }

    }
    //moved left
    else {
            //moved up?
            if (yEnd == (yStart+1)) {

                    double deltaX1 = -(0.5 + x);
                    //lower local origin
                    if(((deltaY / deltaX) * deltaX1 + y) < 0.5) {

                            double deltaY1 = (deltaY / deltaX) * deltaX1;
                            fourBoundaryMove(xStart, yStart, x, y, deltaX1, deltaY1, p, tstep, JxJy);

                            double deltaY2 = 0.5 - y - deltaY1;
                            double deltaX2 = (deltaX1 / deltaY1) * deltaY2;
                            y += deltaY1;
                            fourBoundaryMove(xStart-1, yStart, 0.5, y, deltaX2, deltaY2, p, tstep, JxJy);

                            deltaX -= (deltaX1 + deltaX2);
                            deltaY -= (deltaY1 + deltaY2);
                            x = 0.5 + deltaX2;
                            fourBoundaryMove(xEnd, yEnd, x, -0.5, deltaX, deltaY, p, tstep, JxJy);

                    }
                    //upper local origin
                    else {

                            double deltaY1 = 0.5 - y;
                            deltaX1 = (deltaX / deltaY) * deltaY1;
                            fourBoundaryMove(xStart, yStart, x, y, deltaX1, deltaY1, p, tstep, JxJy);

                            double deltaX2 = -(0.5 + x + deltaX1);
                            double deltaY2 = (deltaY1 / deltaX1) * deltaX2;
                            x += deltaX1;
                            fourBoundaryMove(xStart, yStart+1, x, -0.5, deltaX2, deltaY2, p, tstep, JxJy);

                            deltaX -= (deltaX1 + deltaX2);
                            deltaY -= (deltaY1 + deltaY2);
                            y = deltaY2 - 0.5;
                            fourBoundaryMove(xEnd, yEnd, 0.5, y, deltaX, deltaY,p, tstep, JxJy);
                  
                    }
            }
            //moved down
            else {

                    double deltaY1 = -(0.5 + y);
                    //lower local origin
                    if((-(deltaX / deltaY) * deltaY1 - x) < 0.5) {

                            double deltaX1 = (deltaX / deltaY) * deltaY1;
                            fourBoundaryMove(xStart, yStart, x, y, deltaX1, deltaY1,p, tstep, JxJy);

                            double deltaX2 = -(0.5 + x + deltaX1);
                            double deltaY2 = (deltaY / deltaX) * deltaX2;
                            x += deltaX1;
                            fourBoundaryMove(xStart, yStart-1, x, 0.5, deltaX2, deltaY2, p, tstep, JxJy);

                            deltaX -= (deltaX1 + deltaX2);
                            deltaY -= (deltaY1 + deltaY2);
                            y = 0.5 + deltaY2;
                            fourBoundaryMove(xEnd, yEnd, 0.5, y, deltaX, deltaY, p, tstep, JxJy);

                    }
                    //upper local origin
                    else {

                            double deltaX1 = -(0.5 + x);
                            deltaY1 = (deltaY / deltaX) * deltaX1;
                            fourBoundaryMove(xStart, yStart, x, y, deltaX1, deltaY1, p, tstep, JxJy);

                            double deltaY2 = -(0.5 + y + deltaY1);
                            double deltaX2 = (deltaX1 / deltaY1) * deltaY2;
                            y += deltaY1;
                            fourBoundaryMove(xStart+1, yStart, 0.5, y, deltaX2, deltaY2, p, tstep, JxJy);

                            deltaX -= (deltaX1 + deltaX2);
                            deltaY -= (deltaY1 + deltaY2);
                            x = 0.5 + deltaX2;
                            fourBoundaryMove(xEnd, yEnd, x, 0.5, deltaX, deltaY, p, tstep, JxJy);
                           
                    }
            }
    }
}


__kernel void run_simulation(__global double* particles,
                             __global const double* force, 
                             __global double* cells,
                             __global double* outParticles,
                                  double timeStep,
                                  int n,
                                  double width,
                                  double height,
                                  int numCellsX,//total!!!!!!!!!!!!!!!
                                  int numCellsY,
                                  double cellWidth,
                                  double cellHeight) 
{
   int i = get_global_id(0);
   if((i >= n) || (i % P_SIZE != 0))
        return;


   /*--------------------------------------------------------------------------/
   /---------------- Step 1: particlePush()------------------------------------/
   /--------------------------------------------------------------------------*/

   //a) particle.storePosition() 
        particles[i + PrevX] = particles[i + X];
        particles[i + PrevY] = particles[i + Y];
     
   //b)solver.step(particle, force, timeStep)/----Boris solver-----/
        int j = (int)(i/P_SIZE);

        double getPositionComponentofForceX        = force[j + PositionComponentofForceX];
        double getPositionComponentofForceY        = force[j + PositionComponentofForceY];
        double getBz                               = force[j + Bz];
        double getTangentVelocityComponentOfForceX = force[j + TangentVelocityComponentOfForceX];
        double getTangentVelocityComponentOfForceY = force[j + TangentVelocityComponentOfForceY];
        double getMass                             = particles[i + Mass];

        particles[i + PrevPositionComponentForceX]          = getPositionComponentofForceX;
        particles[i + PrevPositionComponentForceY]          = getPositionComponentofForceY;
        particles[i + PrevBz]                               = getBz;
        particles[i + PrevTangentVelocityComponentOfForceX] = getTangentVelocityComponentOfForceX;
        particles[i + PrevTangentVelocityComponentOfForceY] = getTangentVelocityComponentOfForceY;
   
        double vxminus = particles[i + Vx] + getPositionComponentofForceX * timeStep / (2.0 * getMass);
        double vyminus = particles[i + Vy] + getPositionComponentofForceY * timeStep / (2.0 * getMass);
        
        double t_z = particles[i + Charge] * getBz * timeStep / (2.0 * getMass);
        double s_z = 2 * t_z / (1 + t_z * t_z);

        double vxprime = vxminus + vyminus * t_z;
        double vyprime = vyminus - vxminus * t_z;

        double vxplus = vxminus + vyprime * s_z;
        double vyplus = vyminus - vxprime * s_z;

        particles[i + Vx] = vxplus + getPositionComponentofForceX * timeStep / (2.0 * getMass) + getTangentVelocityComponentOfForceX * timeStep / getMass;
        particles[i + Vy] = vyplus + getPositionComponentofForceY * timeStep / (2.0 * getMass) + getTangentVelocityComponentOfForceY * timeStep / getMass;

        particles[i + X] = particles[i + X] + particles[i + Vx] * timeStep;
        particles[i + Y] = particles[i + Y] + particles[i + Vy] * timeStep;

       /* outParticles[i + Vx] = particles[i + Vx];
        outParticles[i + Vy] = particles[i + Vy];
        outParticles[i + X]  = particles[i + X];
        outParticles[i + Y]  = particles[i + Y];*/

   //c) boundaries.applyOnParticleCenter(solver, force, particle, timeStep)
        double regionBoundaryMapX[9];
        double regionBoundaryMapY[9];
        
        regionBoundaryMapX[X_MIN + Y_MIN] = -width;
        regionBoundaryMapY[X_MIN + Y_MIN] = -height;
        
        regionBoundaryMapX[X_CENTER + Y_MIN] = 0;
        regionBoundaryMapY[X_CENTER + Y_MIN] = -height;

        regionBoundaryMapX[X_MAX + Y_MIN] = width;
        regionBoundaryMapY[X_MAX + Y_MIN] = -height;

        regionBoundaryMapX[X_MIN + Y_CENTER] = -width;
        regionBoundaryMapY[X_MIN + Y_CENTER] = 0;

        regionBoundaryMapX[X_CENTER + Y_CENTER] = 0;
        regionBoundaryMapY[X_CENTER + Y_CENTER] = 0;

        regionBoundaryMapX[X_MAX + Y_CENTER] = width;
        regionBoundaryMapY[X_MAX + Y_CENTER] = 0;

        regionBoundaryMapX[X_MIN + Y_MAX] = -width;
        regionBoundaryMapY[X_MIN + Y_MAX] = height;

        regionBoundaryMapX[X_CENTER + Y_MAX] = 0;
        regionBoundaryMapY[X_CENTER + Y_MAX] = height;

        regionBoundaryMapX[X_MAX + Y_MAX] = width;
        regionBoundaryMapY[X_MAX + Y_MAX] = height;
        
        int regi; 
        regi = get_region(particles[i + X], particles[i + X], particles[i + Y], particles[i + Y], width, height);

        particles[i + X]     = particles[i + X] - regionBoundaryMapX[regi];
        particles[i + PrevX] = particles[i + PrevX] - regionBoundaryMapX[regi];

        particles[i + Y]     = particles[i + Y] - regionBoundaryMapY[regi];
        particles[i + PrevY] = particles[i + PrevY] - regionBoundaryMapY[regi];

        outParticles[i + X] = particles[i + X];

       
   /*--------------------------------------------------------------------------/
   /---------------- Step 2: detector.run()------------------------------------/
   /--------------------------------------------------------------------------*/
   /**Turned off
        int pairsSize = (int)(n/P_SIZE);
        double Pair_a[pairsSize];     //the indexes of the overlapped particles
        double Pair_b[pairsSize];
        int h, k, p = 0;

        for(h = 0; h < n; h += P_SIZE){
            for(k = 0; k < n; k += P_SIZE){
                if(abs(particles[h + X] - particles[k + X]) <= (particles[h + Radius] + particles[k + radius])){
                    if(abs(particles[h + Y] - particles[k + Y]) <= (particles[h + Radius] + particles[k + radius])){
                        Pair_a[p]   = h;
                        Pair_b[p++] = k;
                    }
                }
            }
        }
   */
   /*--------------------------------------------------------------------------/
   /--- Step 3: collisionalgorithm.collide(detector.getOverlappedPairs().....)-/
   /--------------------------------------------------------------------------*/
   /**Turned off-> will get back to this!!!
        int p1, p2;
        double distanceSquare;
        for(k = 0; k < p; k++){
            p1 = Pair_a[k];
            p2 = Pair_b[k];
            
            distanceSquare = ((particles[p1 + X] - particles[p2 + X]) * (particles[p1 + X] - particles[p2 + X]) +
                              (particles[p1 + Y] - particles[p2 + Y]) * (particles[p1 + Y] - particles[p2 + Y]));

            if(distanceSquare <= ((particles[p1 + Radius] + particles[p2 + Radius]) * (particles[p1 + Radius] + particles[p2 + Radius]))){
                s.complete(p1, f, step);
                s.complete(p2, f, step);
                doCollision(p1, p2);
                s.prepare(p1, f, step);
                s.prepare(p2, f, step);
            }
        }
    */

   /*--------------------------------------------------------------------------/
   /---Step 4: interpolation.interpolateToGrid(particles, grid, tstep)---------/
   /--------------------------------------------------------------------------*/
        double JxJy[2]; JxJy[0]=0; JxJy[1]=0; 
        /**X index of local origin i.e. nearest grid point BEFORE particle push*/
        int xStart;
        /**Y index of local origin i.e. nearest grid point BEFORE particle push*/
        int yStart;
        /**X index of local origin i.e. nearest grid point AFTER particle push*/
        int xEnd;
        /**Y index of local origin i.e. nearest grid point AFTER particle push*/
        int yEnd;
        /**Normalized local x coordinate BEFORE particle push*/
        double x;
        /**Normalized local y coordinate BEFORE particle push*/
        double y;
        /**Normalized distance covered in X direction*/
        double deltaX;
        /**Normalized distance covered in X direction*/
        double deltaY;

        //x = particles[i + 
        /*x = p.getPrevX() / g.getCellWidth();
        y = p.getPrevY() / g.getCellHeight();

        xStart = (int) Math.floor(x + 0.5);
        yStart = (int) Math.floor(y + 0.5);

        deltaX = p.getX() / g.getCellWidth();
        deltaY = p.getY() / g.getCellHeight();

        xEnd = (int) Math.floor(deltaX + 0.5);
        yEnd = (int) Math.floor(deltaY + 0.5);

        deltaX -= x;
        deltaY -= y;

        x -= xStart;
        y -= yStart;

      
        //4-boundary move?
        if (xStart == xEnd && yStart == yEnd) {
                fourBoundaryMove(xStart, yStart, x, y, deltaX, deltaY, p, g, tstep);
                }
        //7-boundary move?
        else if (xStart == xEnd || yStart == yEnd) {
                        sevenBoundaryMove(x, y, xStart, yStart, xEnd, yEnd, deltaX, deltaY, p, g, tstep);
                }
                // 10-boundary move
                        else {
                                tenBoundaryMove(x, y, xStart, yStart, xEnd, yEnd, deltaX, deltaY, p, g, tstep);
                        }
                        */
}
