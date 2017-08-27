#ifndef __INTERPOLATE__
#define __INTERPOLATE__

template <typename T>
inline T interpolate(const std::vector<T>& y, T dx, T x, T x0 = (T)0)
{
	int i = (int)floor((x - x0) / dx);
	int j = i + 1;

	if (i < 0)
		return y.front();

	if (j >= (int)y.size())
		return y.back();

	T x_i = x0 + i * dx;
	T x_j = x0 + j * dx;

	T y_i = y[i];
	T y_j = y[j];

	return (x_j*y_i - x_i*y_j + (y_j - y_i)*x) / dx;
}

#endif // __INTERPOLATE__
