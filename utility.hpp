#ifndef __UTILITY_HPP
#define __UTILITY_HPP

namespace bear
{
	template<int SFT, typename T>
	inline T round_shift(T v)
	{
		v = v + (1 << (SFT - 1));

		return v >> (SFT);
	}

	inline int limiteU16(int lum)
	{
		if (lum & 0xffff0000)
		{
			if (lum < 0)return 0;
			return 65535;
		}
		return lum;
	}


	inline int limiteU8(int lum)
	{
		if (lum & 0xffffff00)
		{
			if (lum < 0)return 0;
			return 255;
		}
		return lum;
	}
}

#endif