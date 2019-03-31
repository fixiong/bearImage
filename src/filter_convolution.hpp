#ifndef __FILTERS_CONVOLUTION_HPP
#define __FILTERS_CONVOLUTION_HPP

#include "../../bear/include/image.h"

#include <vector>
#include <utility>
#include <array>

using namespace bear;


template<int Size>
inline size_t size_down(size_t os)
{
	return (os + Size + 1)/2;
}

template<int Size>
inline image_size size_down(image_size sz)
{
	return image_size(size_down<Size>(sz.width),size_down<Size>(sz.height));
}

struct Dfc
{
	void operator()(unsigned short &dst, unsigned int src)
	{
		dst = src;
	}
	void operator()(float &dst, float src)
	{
		dst = src;
	}
};

template<int Size,typename K>
struct KernelADP{};

template<typename K>
struct KernelADP<5,K>
{
	static float run(float v1,float v2,float v3,float v4,float v5)
	{
		return (v1 + v5)*K::w1() + (v2 + v4) * K::w2() + v3 * K::w3();
	}

	static float run_hz(array_ptr<float> p)
	{
		return run(p[0],p[1],p[2],p[3],p[4]);
	}
	static float run_vt(array_ptr<float> p[],size_t i)
	{
		return run(p[0][i],p[1][i],p[2][i],p[3][i],p[4][i]);
	}

	static unsigned int run(
		unsigned short v1,
		unsigned short v2,
		unsigned short v3,
		unsigned short v4,
		unsigned short v5)
	{
		return ((v1 + v5)*K::wi1 + (v2 + v4) * K::wi2 + v3 * K::wi3 + K::ofs) >> K::sft;
	}

	static unsigned int run_hz(array_ptr<unsigned short> p)
	{
		return run(p[0],p[1],p[2],p[3],p[4]);
	}
	static unsigned int run_vt(array_ptr<unsigned short> p[], size_t i)
	{
		return run(p[0][i],p[1][i],p[2][i],p[3][i],p[4][i]);
	}
};


template<typename K>
struct KernelADP<6, K>
{
	static float run(float v1, float v2, float v3, float v4, float v5, float v6)
	{
		return (v1 + v6)*K::w1() + (v2 + v5) * K::w2() + (v3 + v4) * K::w3();
	}

	static float run_hz(array_ptr<float> p)
	{
		return run(p[0], p[1], p[2], p[3], p[4], p[5]);
	}
	static float run_vt(array_ptr<float> p[], size_t i)
	{
		return run(p[0][i], p[1][i], p[2][i], p[3][i], p[4][i], p[5][i]);
	}

	static unsigned int run(
		unsigned short v1,
		unsigned short v2,
		unsigned short v3,
		unsigned short v4,
		unsigned short v5,
		unsigned short v6)
	{
		return ((v1 + v6)*K::wi1 + (v2 + v5) * K::wi2 + (v3 + v4) * K::wi3 + K::ofs) >> K::sft;
	}

	static unsigned int run_hz(array_ptr<unsigned short> p)
	{
		return run(p[0], p[1], p[2], p[3], p[4], p[5]);
	}
	static unsigned int run_vt(array_ptr<unsigned short> p[], size_t i)
	{
		return run(p[0][i], p[1][i], p[2][i], p[3][i], p[4][i], p[5][1]);
	}
};

/////////////////////////////////////////////////////////////////////

template<int Size, typename Kernel>
struct _UpFilter
{
};

template<typename K, typename C>
inline void upf5(
	array_ptr<float> srow1,
	array_ptr<float> srow2,
	array_ptr<float> srow3,
	array_ptr<float> drow,
	array_ptr<float> drow1,
	size_t &x,
	C &&c)
{
	float s12_32 = srow1[2] + srow3[2];
	float s11_31 = srow1[1] + srow3[1];
	float s21_s31 = srow2[1] + srow3[1];
	float s22_s32 = srow2[2] + srow3[2];
	std::forward<C>(c)(drow[x], (srow1[0] + srow3[0] + s12_32) * (K::w1()*K::w1()) +
		(srow2[0] + srow2[2] + s11_31) * (K::w1()*K::w3()) +
		srow2[1] * (K::w3()*K::w3()));
	std::forward<C>(c)(drow1[x], (srow2[0] + srow3[0] + s22_s32) * (K::w1()*K::w2()) +
		(s21_s31) * (K::w2()*K::w3()));
	++x;
	if (x >= (int)drow.size())return;
	std::forward<C>(c)(drow[x], (s11_31 + s12_32) * (K::w1()*K::w2()) +
		(srow2[1] + srow2[2]) * (K::w2()*K::w3()));
	std::forward<C>(c)(drow1[x], (s21_s31 + s22_s32) * (K::w2()*K::w2()));
}


template<typename K, typename C>
inline void upf5(
	array_ptr<unsigned short> srow1,
	array_ptr<unsigned short> srow2,
	array_ptr<unsigned short> srow3,
	array_ptr<unsigned short> drow,
	array_ptr<unsigned short> drow1,
	size_t &x,
	C &&c)
{
	unsigned int s12_32 = srow1[2] + srow3[2];
	unsigned int s11_31 = srow1[1] + srow3[1];
	unsigned int s21_s31 = srow2[1] + srow3[1];
	unsigned int s22_s32 = srow2[2] + srow3[2];
	std::forward<C>(c)(drow[x], ((srow1[0] + srow3[0] + s12_32) * ((K::wi1*K::wi1 + (1 << (K::sft - 1))) >> K::sft) +
		(srow2[0] + srow2[2] + s11_31) * ((K::wi1*K::wi3 + (1 << (K::sft - 1))) >> K::sft) +
		srow2[1] * ((K::wi3*K::wi3 + (1 << (K::sft - 1))) >> K::sft) + K::ofs) >> K::sft);
	std::forward<C>(c)(drow1[x], ((srow2[0] + srow3[0] + s22_s32) * ((K::wi1*K::wi2 + (1 << (K::sft - 1))) >> K::sft) +
		(s21_s31) * ((K::wi2*K::wi3 + (1 << (K::sft - 1))) >> K::sft) + K::ofs) >> K::sft);
	++x;
	if (x >= (int)drow.size())return;
	std::forward<C>(c)(drow[x], ((s11_31 + s12_32) * ((K::wi1*K::wi2 + (1 << (K::sft - 1))) >> K::sft) +
		(srow2[1] + srow2[2]) * ((K::wi2*K::wi3 + (1 << (K::sft - 1))) >> K::sft) + K::ofs) >> K::sft);
	std::forward<C>(c)(drow1[x], ((s21_s31 + s22_s32) * ((K::wi2*K::wi2 + (1 << (K::sft - 1))) >> K::sft) + K::ofs) >> K::sft);
}


template<typename K>
struct _UpFilter<5, K>
{
	template<typename Image, typename C = Dfc>
	static void run(Image dst, Image src, C &&c = Dfc())
	{
		using Unit = typename Image::elm_type;

		array_ptr<Unit> srow1, srow2, srow3;
		size_t dw = dst.width();
		array_ptr<Unit> drow1;

		std::vector<Unit> ct;
		for (size_t y = 0; y<dst.height(); ++y)
		{
			size_t sy = y >> 1;
			srow1 = src[sy];
			srow2 = src[sy + 1];
			srow3 = src[sy + 2];
			array_ptr<Unit> drow = dst[y];

			++y;
			if (y >= dst.height())
			{
				ct.resize(dw);
				drow1 = ct;
			}
			else drow1 = dst[y];


			for (size_t x = 0; x<dw; ++x)
			{
				upf5<K>(
					srow1.clip(x,srow1.size()),
					srow2.clip(x,srow2.size()),
					srow3.clip(x,srow3.size()),
					drow, drow1, x, std::forward<C>(c));
			}
		}
	}
};


template<typename Kernel>
struct UpFilter :public _UpFilter<Kernel::Size, Kernel> {};


////////////////////////////////////////////////////////////////////////////

template<int Size,typename Kernel>
struct DDKernel
{};

template<typename Kernel>
struct DDKernel<5, Kernel>
{
	enum { Size = 6 };
	enum { 
		wi1 = Kernel::wi1,
		wi2 = Kernel::wi1 + Kernel::wi2,
		wi3 = Kernel::wi2 + Kernel::wi3,
		ofs = Kernel::ofs, sft = Kernel::sft };

	constexpr static float w1()
	{
		return Kernel::w1();
	}
	constexpr static float w2()
	{
		return Kernel::w1() + Kernel::w2();
	}
	constexpr static float w3()
	{
		return Kernel::w2() + Kernel::w3();
	}
};

template<typename T>
struct _ZERO_UNIT{};

template<>
struct _ZERO_UNIT<float>
{
	constexpr static float run()
	{
		return 0.0f;
	}

	static float to_zero(float v)
	{
		return v;
	}
};

template<>
struct _ZERO_UNIT<unsigned short>
{
	constexpr static unsigned short run()
	{
		return 32768;
	}


	static unsigned short to_zero(unsigned short v)
	{
		return v - 32768;
	}
};


template<typename _K, typename Line>
void down_line_dx(Line dst,Line src)
{
	using Unit = typename Line::elm_type;

	typedef DDKernel<_K::Size,_K> K;
	typedef _ZERO_UNIT<Unit> zu;

	auto sw = src.size();
	auto dw = size_down<_K::Size>(sw);

	std::array<Unit,K::Size> tmp;

	size_t lm = K::Size >> 1;

	int s = 1 - K::Size;

	for (size_t i = 0; i < lm; ++i)
	{
		for (int j = 0; j<K::Size; ++j)
		{
			if (j + s<0)
			{
				tmp[j] = _ZERO_UNIT<Unit>::run();
			}
			else
			{
				tmp[j] = src[j + s];
			}
		}
		dst[i] = zu::to_zero(KernelADP<K::Size, K>::run_hz(array_ptr<Unit>(tmp)));
		s += 2;
	}

	auto mi = sw >> 1;

	for (auto i = lm; i<mi; ++i)
	{
		dst[i] = zu::to_zero(KernelADP<K::Size, K>::run_hz(src.clip(s,src.size())));
		s += 2;
	}

	for (auto i = mi; i<dw; ++i)
	{
		for (size_t j = 0; j<K::Size; ++j)
		{
			if (j + s >= sw)
			{
				tmp[j] = _ZERO_UNIT<Unit>::run();
			}
			else
			{
				tmp[j] = src[j + s];
			}
		}
		dst[i] = zu::to_zero(KernelADP<K::Size, K>::run_hz(array_ptr<Unit>(tmp)));
		s += 2;
	}
}

template<typename K, typename Line>
void down_line(Line dst, Line src)
{
	using Unit = typename Line::elm_type;

	auto sw = src.size();
	auto dw = size_down<K::Size>(sw);

	std::array<Unit,K::Size> tmp;

	size_t lm = K::Size >> 1;

	int s = 1 - K::Size;

	for (size_t i = 0; i < lm; ++i)
	{
		for (size_t j = 0; j<K::Size; ++j)
		{
			if (j + s<0)
			{
				tmp[j] = src[0];
			}
			else
			{
				tmp[j] = src[j + s];
			}
		}
		dst[i] = KernelADP<K::Size, K>::run_hz(array_ptr<Unit>(tmp));
		s += 2;
	}

	auto mi = sw >> 1;

	for (auto i = lm; i<mi; ++i)
	{
		dst[i] = KernelADP<K::Size, K>::run_hz(src.clip(s,src.size()));
		s += 2;
	}

	for (auto i = mi; i<dw; ++i)
	{
		for (size_t j = 0; j<K::Size; ++j)
		{
			if (j + s >= sw)
			{
				tmp[j] = src[sw - 1];
			}
			else
			{
				tmp[j] = src[j + s];
			}
		}
		dst[i] = KernelADP<K::Size, K>::run_hz(array_ptr<Unit>(tmp));
		s += 2;
	}
}

template<int Size, typename Line>
inline void roll_buf(Line buf[])
{
	auto tmp1 = buf[0];
	auto tmp2 = buf[1];

	for (size_t i = 0; i<Size - 2; ++i)
	{
		buf[i] = buf[i + 2];
	}

	buf[Size - 2] = tmp1;
	buf[Size - 1] = tmp2;
}

template<bool D,typename K>
struct _DOWN_LINE
{
	template<typename Line>
	static void run(Line dst, Line src)
	{
		down_line<K>(dst, src);
	}
};


template<typename K>
struct _DOWN_LINE<true,K>
{
	template<typename Line>
	static void run(Line dst, Line src)
	{
		down_line_dx<K>(dst, src);
	}
};

template<int Size, typename K>
struct _DownFilter
{
	template<bool D, typename Image, typename C = Dfc>
	static void _run(Image dst, Image src, C &&c = Dfc())
	{
		using Unit = typename Image::elm_type;

		typedef _DOWN_LINE<D,K> DL;


		image<Unit,1> bufs(dst.width(), Size);
		array_ptr<Unit> buf[Size];
		for (int i = 0; i<Size; ++i)
		{
			buf[i] = bufs[i];
		}

		for (int j = 0; j<Size - 2; ++j)
		{
			DL::run(buf[j], src[0]);
		}

		size_t sh = src.height();

		size_t s = - 1;

		for (size_t i = 0; i<dst.height(); ++i)
		{


			for (size_t j = Size - 2; j<Size; ++j)
			{
				if (s < 0)
				{
					DL::run(buf[j], src[0]);
				}
				else if (s >= sh)
				{
					DL::run(buf[j], src[sh - 1]);
				}
				else
				{
					DL::run(buf[j], src[s]);
				}
				++s;
			}

			size_t dw = dst.width();
			auto drow = dst[i];

			for (size_t j = 0; j<dw; ++j)
			{
				std::forward<C>(c)(drow[j], KernelADP<Size, K>::run_vt(buf, j));
			}

			roll_buf<Size>(buf);
		}
	}

	template<typename Image, typename C = Dfc>
	static void run(Image dst, Image src, C &&c = Dfc())
	{
		_run<false>(dst, src, c);
	}

	template<typename Image, typename C = Dfc>
	static void run_dx(Image dst, Image src, C &&c = Dfc())
	{
		_run<true>(dst, src, c);
	}

	template<typename Image, typename C = Dfc>
	static void run_dy(Image dst, Image src, C &&c = Dfc())
	{
		using Unit = typename Image::elm_type;

		typedef _DOWN_LINE<false, K> DL;

		typedef DDKernel<K::Size, K> _K;
		const int _Size = _K::Size;

		image<Unit,1> bufs(dst.width(), _Size);
		array_ptr<Unit> buf[_Size];
		for (int i = 0; i<_Size; ++i)
		{
			buf[i] = bufs[i];
		}

		for (int j = 0; j<_Size - 2; ++j)
		{
			buf[j].fill(_ZERO_UNIT<Unit>::run());
		}

		auto sh = src.height();

		size_t s = -1;

		for (size_t i = 0; i<dst.height(); ++i)
		{


			for (size_t j = _Size - 2; j<_Size; ++j)
			{
				if (s < 0)
				{
					buf[j].fill(_ZERO_UNIT<Unit>::run());
				}
				else if (s >= sh)
				{
					buf[j].fill(_ZERO_UNIT<Unit>::run());
				}
				else
				{
					DL::run(buf[j], src[s]);
				}
				++s;
			}

			size_t dw = dst.width();
			auto drow = dst[i];

			for (size_t j = 0; j<dw; ++j)
			{
				std::forward<C>(c)(drow[j], _ZERO_UNIT<Unit>::to_zero(KernelADP<_Size, _K>::run_vt(buf, j)));
			}

			roll_buf<_Size>(buf);
		}
	}
};

template<typename Kernel>
struct DownFilter :public _DownFilter<Kernel::Size, Kernel> {};




#endif
