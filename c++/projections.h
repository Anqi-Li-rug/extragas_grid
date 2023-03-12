/* projections.h
   .
   PROJECTIONS IS A CLASS OF TRASFORMATIONS THAT REDUCE THE NUMBER OF 
   COORDINATES OF THE PHASE SPACE OF AN AXI-SYMMETRIC SYSTEM.
   .
   TWO TRASFORMATIONS:
   - EXTERNAL VIEW: 
   The observer is outside the system, the line of nodes is along Y, 
   the inclination (incl) of the symmetry axis and the line of sight is a 
   free parameter.
   - INTERNAL VIEW:
   The observer is inside the system, on the plane of symmetry, at a 
   distance R0 from the centre and moving at a velocity v0_phi. The 
   azimuthal (l) coordinate is set at 0 towards the centre of the system.
 */

#include "array.h"
#include "findValue.h"


using namespace ffarray;

class projections {
 private:
  double R, phi, z, v_R, v_phi, v_z; // [kpc], [km/s]
  float xi, yi, vrad;
  float l, b, d, vlos, hvlos;
  float x_proj(float, float);
  float y_proj(float, float, float, float);
  float v_proj(float, float, float, float, float);
  float longitude(float, float, float);
  float latitude(float, float, float, float, float);
  float v_los(float, float, float, float, float, float);
  float distance(float, float, float, float, float);
 public:
  projections();
  projections(double, double, double, double, double, double);
  void external_view(float,float, float&, float&, float& );
  void internal_view(float, float, float&, float&, float&, float& );
  void internal_view_hvcs(float, float, float, float, float&, float&, float&, float& );
  float modelgal(float, float, int, float, float, float, float);
  float hvcvel(float, float, float, float, float, float, float);

  void BuildInternalCube(Array3D<float>& , Array2D<float>&,
			 int, int, int, float, float, float);

  ~projections();
};
