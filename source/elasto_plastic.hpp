#ifndef __ELASTO_PLASTIC__
#define __ELASTO_PLASTIC__

//
//       ^ d_y
//       |             x  new displacement
//       |            /
//       |     ------/ 
//       |   /      O  \  new origin
//       |  |      /    |
//       | |      o      |
//       |  |   origin  |
//       |   \         /
//       |     -------      d_x
//    ------------------------>
//       |

inline void elasto_plastic(
		double k,
		double d_yield,
		double& d_o_x,
		double& d_o_y,
		double d_x,
		double d_y,
		double& F_x,
		double& F_y,
		double& k_tan)
{
	double F_yield = k * d_yield;

	double delta_x = d_x - d_o_x;
	double delta_y = d_y - d_o_y;
	double delta_r = sqrt(sq(delta_x) + sq(delta_y));

	if (delta_r < d_yield)  // elastic
	{
		F_x = k * delta_x;
		F_y = k * delta_y;
		
		k_tan = k;
	}
	else  // plastic
	{
		F_x = (delta_x / delta_r) * F_yield;
		F_y = (delta_y / delta_r) * F_yield;

		d_o_x = d_x - (d_yield / delta_r) * delta_x;
		d_o_y = d_y - (d_yield / delta_r) * delta_y;

		k_tan = 0.0;
	}
}


/* 1D model with stiffness degradation
inline void elasto_plastic(
		double k_0,
		double d_yield,
		double c,
		double& d_o,
		double d,
		double& F,
		double& k_tan)
{
	double F_yield = k_0 * d_yield;
	double delta = d - d_o;

	double k = (c - 2.0) * F_yield / ((c - 2.0) * d_yield - c * abs(d_o));

	if (abs(delta) < F_yield / k)  // elastic
	{
		F = k * delta;
		k_tan = k;
	}
	else  // plastic
	{
		if (delta > 0.0)  // positive force
		{
			if (d >= d_yield)  // origin on the right
			{
				double a = d - d_yield;
				k = 2.0 * F_yield / (c * a + 2.0 * d_yield);
				d_o = d - F_yield / k;
				F = F_yield;
			}
			else  // origin on the left
			{
				double a = (-d + d_yield) / (1.0 - c);
				k = 2.0 * F_yield / (c * a + 2.0 * d_yield);
				d_o = d - F_yield / k;
				F = F_yield;
			}
		}
		else  // negative force
		{
			if (d <= -d_yield)  // origin on the left
			{
				double a = -d - d_yield;
				k = 2.0 * F_yield / (c * a + 2.0 * d_yield);
				d_o = d + F_yield / k;
				F = -F_yield;
			}
			else  // origin on the right
			{
				double a = (d + d_yield) / (1.0 - c);
				k = 2.0 * F_yield / (c * a + 2.0 * d_yield);
				d_o = d + F_yield / k;
				F = -F_yield;
		
			}
		}

		k_tan = 0.0;
	}
}
*/

#endif // __ELASTO_PLASTIC__
