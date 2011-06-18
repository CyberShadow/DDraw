module draw;

import std.file;
import std.string;
import std.algorithm;
import std.exception;
import std.conv;
import std.math;
import std.traits;
static import core.bitop;

struct Image(COLOR)
{
	alias COLOR ColorType;

	uint w, h;
	COLOR[] pixels;

	this(uint w, uint h)
	{
		size(w, h);
	}

	void size(uint w, uint h)
	{
		if ((this.w && this.h) && (w && h))
			throw new Exception("Resize not implemented");
		this.w = w;
		this.h = h;
		pixels.length = w*h;
	}

	COLOR opIndex(uint x, uint y)
	{
		return pixels[y*w+x];
	}

	void opIndexAssign(COLOR value, uint x, uint y)
	{
		pixels[y*w+x] = value;
	}

	void hline(uint x1, uint x2, uint y, COLOR c)
	{
		auto rowOffset = y*w;
		pixels[rowOffset+x1..rowOffset+x2] = c;
	}

	void vline(uint x, uint y1, uint y2, COLOR c)
	{
		foreach (y; y1..y2) // TODO: optimize
			pixels[y*w+x] = c;
	}

	void rect(uint x1, uint y1, uint x2, uint y2, COLOR c) // []
	{
		hline(x1, x2+1, y1, c);
		hline(x1, x2+1, y2, c);
		vline(x1, y1, y2+1, c);
		vline(x2, y1, y2+1, c);
	}

	void fillRect(uint x1, uint y1, uint x2, uint y2, COLOR b) // [)
	{
		foreach (y; y1..y2)
			pixels[y*w+x1..y*w+x2] = b;
	}

	void fillRect(uint x1, uint y1, uint x2, uint y2, COLOR c, COLOR b) // []
	{
		rect(x1, y1, x2, y2, c);
		foreach (y; y1+1..y2)
			pixels[y*w+x1+1..y*w+x2] = b;
	}

	static string[] readPNMHeader(ref ubyte[] data)
	{
		string[] fields;
		uint wordStart = 0;
		uint p;
		for (p=1; p<data.length && fields.length<4; p++)
			if (!iswhite(data[p-1]) && iswhite(data[p]))
				fields ~= cast(string)data[wordStart..p];
			else
			if (iswhite(data[p-1]) && !iswhite(data[p]))
				wordStart = p;
		data = data[p..$];
		enforce(fields.length==4, "Header too short");
		enforce(fields[0].length==2 && fields[0][0]=='P', "Invalid signature");
		return fields;
	}

	void savePNM()(string filename) // RGB only
	{
		static assert(is(typeof(COLOR.init.r)) && is(typeof(COLOR.init.g)) && is(typeof(COLOR.init.b)) && !is(typeof(COLOR.init.a)), "PNM only supports RGB");
		static assert(COLOR.init.r.sizeof == COLOR.init.g.sizeof && COLOR.init.g.sizeof == COLOR.init.b.sizeof, "Inconsistent color channel sizes");
		alias typeof(COLOR.init.r) CHANNEL_TYPE;
		enforce(w*h == pixels.length, "Dimension mismatch");
		ubyte[] header = cast(ubyte[])format("P6\n%d %d %d\n", w, h, CHANNEL_TYPE.max);
		ubyte[] data = new ubyte[header.length + pixels.length * 3 * CHANNEL_TYPE.sizeof];
		data[0..header.length] = header;
		CHANNEL_TYPE* p = cast(CHANNEL_TYPE*)(data.ptr + header.length);
		foreach (c; pixels)
		{
			*p++ = bswap(c.r);
			*p++ = bswap(c.g);
			*p++ = bswap(c.b);
		}
		std.file.write(filename, data);
	}

	void loadPGM()(string filename)
	{
		static assert(is(typeof(COLOR.init == 0)), "PGM only supports grayscale");
		ubyte[] data = cast(ubyte[])read(filename);
		string[] fields = readPNMHeader(data);
		enforce(fields[0]=="P5", "Invalid signature");
		w = to!uint(fields[1]);
		h = to!uint(fields[2]);
		enforce(data.length / COLOR.sizeof == w*h, "Dimension / filesize mismatch");
		enforce(to!uint(fields[3]) == COLOR.max);
		pixels = cast(COLOR[])data;
		static if (COLOR.sizeof > 1)
			foreach (ref pixel; pixels)
				pixel = bswap(pixel);
	}

	void savePGM()(string filename)
	{
		static assert(is(typeof(COLOR.init == 0)), "PGM only supports grayscale");
		ubyte[] header = cast(ubyte[])format("P5\n%d %d\n%d\n", w, h, COLOR.max);
		ubyte[] data = new ubyte[header.length + pixels.length * COLOR.sizeof];
		data[0..header.length] = header;
		COLOR* p = cast(COLOR*)(data.ptr + header.length);
		foreach (c; pixels)
			*p++ = bswap(c);
		std.file.write(filename, data);
	}

	void loadRGBA()(string filename, uint w, uint h)
	{
		static assert(is(typeof(COLOR.init.r)) && is(typeof(COLOR.init.g)) && is(typeof(COLOR.init.b)) && is(typeof(COLOR.init.a)), "COLOR is not RGBA");
		pixels = cast(COLOR[])read(filename);
		enforce(pixels.length == w*h, "Dimension / filesize mismatch");
		this.w = w;
		this.h = h;
	}

	static T bswap(T)(T b)
	{
		static if (b.sizeof == 1)
			return b;
		else
		static if (b.sizeof == 2)
			return cast(T)((b >> 8) | (b << 8));
		else
		static if (b.sizeof == 4)
			return core.bitop.bswap(b);
		else
			static assert(false, "Don't know how to bswap " ~ T.stringof);
	}

	Image!COLOR crop(uint x1, uint y1, uint x2, uint y2)
	{
		auto nw = x2-x1;
		auto nh = y2-y1;
		auto newImage = Image!COLOR(nw, nh);
		uint oldOffset, newOffset;
		foreach (y; y1..y2)
		{
			auto oldOffset2 = oldOffset + w;
			auto newOffset2 = newOffset + nw;
			newImage.pixels[newOffset..newOffset2] = pixels[oldOffset..oldOffset2];
			oldOffset = oldOffset2;
			newOffset = newOffset2;
		}
		return newImage;
	}

	Image!COLOR2 convert(string pred=`c`, COLOR2=typeof(((){ COLOR c; size_t i; return mixin(pred); })()))()
	{
		auto newImage = Image!COLOR2(w, h);
		foreach (i, c; pixels)
			newImage.pixels[i] = mixin(pred);
		return newImage;
	}
}

struct RGB  { ubyte r, g, b, x; }
struct RGBA { ubyte r, g, b, a; }
struct GA   { ubyte g, a; }
struct GA16 { ushort g, a; }

// *****************************************************************************

struct HRImage(COLOR, uint HR)
{
	Image!COLOR hr, lr;

	this(uint w, uint h)
	{
		lr.size(w, h);
		hr.size(w*HR, h*HR);
	}

	void upscale()
	{
		foreach (y; 0..lr.h)
			foreach (x, c; lr.pixels[y*lr.w..(y+1)*lr.w])
				hr.fillRect(x*HR, y*HR, x*HR+HR, y*HR+HR, c);
	}

	void downscale()
	{
		foreach (y; 0..lr.h)
			foreach (x; 0..lr.w)
			{
				static assert(HR*HR <= 256);
				static if (is(typeof(COLOR.init.a))) // downscale with alpha
				{
					ExpandType!(COLOR, 1+COLOR.init.a.sizeof) sum;
					ExpandType!(typeof(COLOR.init.a), 1) alphaSum;
					auto start = y*HR*hr.w + x*HR;
					foreach (j; 0..HR)
					{
						foreach (p; hr.pixels[start..start+HR])
						{
							foreach (i, f; p.tupleof)
								static if (p.tupleof[i].stringof != "p.a")
								{
									enum FIELD = p.tupleof[i].stringof[2..$];
									mixin("sum."~FIELD~" += cast(typeof(sum."~FIELD~"))p."~FIELD~" * p.a;");
								}
							alphaSum += p.a;
						}
						start += hr.w;
					}
					if (alphaSum)
					{
						auto result = cast(COLOR)(sum / alphaSum);
						result.a = cast(typeof(result.a))(alphaSum / (HR*HR));
						lr[x, y] = result;
					}
					else
					{
						static assert(COLOR.init.a == 0);
						lr[x, y] = COLOR.init;
					}
				}
				else
				{
					ExpandType!(COLOR, 1) sum;
					auto start = y*HR*hr.w + x*HR;
					foreach (j; 0..HR)
					{
						foreach (p; hr.pixels[start..start+HR])
							sum += p;
						start += hr.w;
					}
					lr[x, y] = cast(COLOR)(sum / (HR*HR));
				}
			}
	}

	void pixel(uint x, uint y, COLOR c)
	{
		pixelHR(x*HR, y*HR, c);
	}
	
	void pixelHR(uint x, uint y, COLOR c)
	{
		hr.fillRect(x, y, x+HR, y+HR, c);
	}

	void line(uint x1, uint y1, uint x2, uint y2, COLOR c)
	{
		auto xmin = min(x1, x2);
		auto xmax = max(x1, x2);
		auto ymin = min(y1, y2);
		auto ymax = max(y1, y2);

		if (xmax-xmin > ymax-ymin)
			foreach (x; min(x1,x2)..max(x1,x2)+1)
				pixelHR(x*HR, itpl(y1*HR, y2*HR, x, x1, x2), c);
		else
			foreach (y; min(y1,y2)..max(y1,y2)+1)
				pixelHR(itpl(x1*HR, x2*HR, y, y1, y2), y*HR, c);
	}
}

template ExpandType(T, uint BYTES)
{
	static if (is(T : ulong))
	{
		static if (T.sizeof + BYTES <= 2)
			alias ushort ExpandType;
		else
		static if (T.sizeof + BYTES <= 4)
			alias uint ExpandType;
		else
		static if (T.sizeof + BYTES <= 8)
			alias ulong ExpandType;
		else
			static assert(0, "No type big enough to fit " ~ T.sizeof.stringof ~ " + " ~ BYTES.stringof ~ " bytes");
	}
	else
	static if (is(T==struct))
		struct ExpandType
		{
			static string mixFields()
			{
				string s;
				string[] fields;
				foreach (i, f; T.init.tupleof)
				{
					string field = T.tupleof[i].stringof;
					while (field[0] != '.')
						field = field[1..$];
					field = field[1..$];
					if (field != "x") // HACK
						fields ~= field;
				}

				foreach (field; fields)
					s ~= "ExpandType!(typeof(" ~ T.stringof ~ ".init." ~ field ~ "), "~BYTES.stringof~") " ~ field ~ ";\n";
				s ~= "\n";

				s ~= "void opOpAssign(string OP)(" ~ T.stringof ~ " color) if (OP==`+`)\n";
				s ~= "{\n";
				foreach (field; fields)
					s ~= "	"~field~" += color."~field~";\n";
				s ~= "}\n\n";

				s ~= T.stringof ~ " opBinary(string OP, T)(T divisor) if (OP==`/`)\n";
				s ~= "{\n";
				s ~= "	"~T.stringof~" color;\n";
				foreach (field; fields)
					s ~= "	color."~field~" = cast(typeof(color."~field~")) ("~field~" / divisor);\n";
				s ~= "	return color;\n";
				s ~= "}\n\n";

				return s;
			}

			//pragma(msg, mixFields());
			mixin(mixFields());
		}
	else
		static assert(0);
}

// *****************************************************************************

T itpl(T, U)(T low, T high, U r, U rLow, U rHigh)
{
	return cast(T)(low + (cast(Signed!T)high-cast(Signed!T)low) * (cast(Signed!U)r - cast(Signed!U)rLow) / (cast(Signed!U)rHigh - cast(Signed!U)rLow));
}

// *****************************************************************************

struct GammaRamp(LUM_COLOR, PIX_COLOR)
{
	LUM_COLOR[PIX_COLOR.max+1] pix2lum;
	PIX_COLOR[LUM_COLOR.max+1] lum2pix;

	this(double gamma)
	{
		foreach (pix; 0..PIX_COLOR.max+1)
			pix2lum[pix] = cast(LUM_COLOR)(pow(pix/cast(double)PIX_COLOR.max,   gamma)*LUM_COLOR.max);
		foreach (lum; 0..LUM_COLOR.max+1)
			lum2pix[lum] = cast(PIX_COLOR)(pow(lum/cast(double)LUM_COLOR.max, 1/gamma)*PIX_COLOR.max);
	}

	Image!PIX_COLOR image2pix(in Image!LUM_COLOR lumImage)
	{
		auto pixImage = Image!PIX_COLOR(lumImage.w, lumImage.h);
		foreach (i, p; lumImage.pixels)
			pixImage.pixels[i] = lum2pix[p];
		return pixImage;
	}
}
