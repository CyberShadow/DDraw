module draw;

import std.file;
import std.string;
import std.algorithm;
import std.exception;
import std.conv;
import std.math;
import std.traits;
import std.zlib;
import crc32;
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

	ref COLOR opIndex(uint x, uint y)
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

	// Unchecked! Make sure area is bounded.
	void uncheckedFloodFill(uint x, uint y, COLOR c)
	{
		floodFillPtr(&this[x, y], c, this[x, y]);
	}

	private void floodFillPtr(COLOR* pp, COLOR c, COLOR f)
	{
		COLOR* p0 = pp; while (*p0==f) p0--; p0++;
		COLOR* p1 = pp; while (*p1==f) p1++; p1--;
		for (auto p=p0; p<=p1; p++)
			*p = c;
		p0 -= w; p1 -= w;
		for (auto p=p0; p<=p1; p++)
			if (*p == f)
				floodFillPtr(p, c, f);
		p0 += w*2; p1 += w*2;
		for (auto p=p0; p<=p1; p++)
			if (*p == f)
				floodFillPtr(p, c, f);
	}

	template SameSize(T, U...)
	{
		static if (U.length)
			enum SameSize = T.sizeof==U[0].sizeof && SameSize!U;
		else
			enum SameSize = true;
	}

	template ChannelType(T)
	{
		static if (is(T : ulong))
			alias T ChannelType;
		else
		static if (is(T == struct))
		{
			static assert(SameSize!T, "Inconsistent color channel sizes");
			alias typeof(T.init.tupleof[0]) ChannelType;
		}
		else
			static assert(0, "Can't get channel type of " ~ T.stringof);
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
		static assert(__traits(allMembers, COLOR).stringof == `tuple("r","g","b")`, "PNM only supports RGB");
		alias ChannelType!COLOR CHANNEL_TYPE;
		enforce(w*h == pixels.length, "Dimension mismatch");
		ubyte[] header = cast(ubyte[])format("P6\n%d %d %d\n", w, h, CHANNEL_TYPE.max);
		ubyte[] data = new ubyte[header.length + pixels.length * COLOR.sizeof];
		data[0..header.length] = header;
		data[header.length..$] = cast(ubyte[])pixels;
		static if (CHANNEL_TYPE.sizeof > 1)
		{
			auto end = cast(CHANNEL_TYPE*)data.ptr+data.length;
			for (CHANNEL_TYPE* p = cast(CHANNEL_TYPE*)(data.ptr + header.length); p<end; p++)
				*p = bswap(*p);
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
		static assert(__traits(allMembers, COLOR).stringof == `tuple("r","g","b","a")`, "COLOR is not RGBA");
		pixels = cast(COLOR[])read(filename);
		enforce(pixels.length == w*h, "Dimension / filesize mismatch");
		this.w = w;
		this.h = h;
	}

	void savePNG()(string filename)
	{
		enum : ulong { SIGNATURE = 0x0a1a0a0d474e5089 }

		struct PNGChunk
		{
			char[4] type;
			const(void)[] data;

			uint crc32()
			{
				uint crc = strcrc32(type);
				foreach (v; cast(ubyte[])data)
					crc = update_crc32(v, crc);
				return ~crc;
			}

			this(string type, const(void)[] data)
			{
				this.type[] = type;
				this.data = data;
			}
		}

		enum PNGColourType : ubyte { G, RGB=2, PLTE, GA, RGBA=6 }
		enum PNGCompressionMethod : ubyte { DEFLATE }
		enum PNGFilterMethod : ubyte { ADAPTIVE }
		enum PNGInterlaceMethod : ubyte { NONE, ADAM7 }

		enum PNGFilterAdaptive : ubyte { NONE, SUB, UP, AVERAGE, PAETH }

		struct PNGHeader
		{
		align(1):
			uint width, height;
			ubyte colourDepth;
			PNGColourType colourType;
			PNGCompressionMethod compressionMethod;
			PNGFilterMethod filterMethod;
			PNGInterlaceMethod interlaceMethod;
			static assert(PNGHeader.sizeof == 13);
		}

		alias ChannelType!COLOR CHANNEL_TYPE;

		static if (!is(COLOR == struct))
			enum COLOUR_TYPE = PNGColourType.G;
		else
		static if (structFields!COLOR()==["r","g","b"])
			enum COLOUR_TYPE = PNGColourType.RGB;
		else
		static if (structFields!COLOR()==["g","a"])
			enum COLOUR_TYPE = PNGColourType.GA;
		else
		static if (structFields!COLOR()==["r","g","b","a"])
			enum COLOUR_TYPE = PNGColourType.RGBA;
		else
			static assert(0, "Unsupported PNG color type: " ~ COLOR.stringof);

		PNGChunk[] chunks;
		PNGHeader header = {
			width : bswap(w),
			height : bswap(h),
			colourDepth : CHANNEL_TYPE.sizeof * 8,
			colourType : COLOUR_TYPE,
			compressionMethod : PNGCompressionMethod.DEFLATE,
			filterMethod : PNGFilterMethod.ADAPTIVE,
			interlaceMethod : PNGInterlaceMethod.NONE,
		};
		chunks ~= PNGChunk("IHDR", cast(void[])[header]);
		uint idatStride = w*COLOR.sizeof+1;
		ubyte[] idatData = new ubyte[h*idatStride];
		for (uint y=0; y<h; y++)
		{
			idatData[y*idatStride] = PNGFilterAdaptive.NONE;
			auto rowPixels = cast(COLOR[])idatData[y*idatStride+1..(y+1)*idatStride];
			rowPixels[] = pixels[y*w..(y+1)*w];

			static if (CHANNEL_TYPE.sizeof > 1)
				foreach (ref p; cast(CHANNEL_TYPE[])rowPixels)
					p = bswap(p);
		}
		chunks ~= PNGChunk("IDAT", compress(idatData, 5));
		chunks ~= PNGChunk("IEND", null);

		uint totalSize = 8;
		foreach (chunk; chunks)
			totalSize += 8 + chunk.data.length + 4;
		ubyte[] data = new ubyte[totalSize];

		*cast(ulong*)data.ptr = SIGNATURE;
		uint pos = 8;
		foreach(chunk;chunks)
		{
			uint i = pos;
			uint chunkLength = chunk.data.length;
			pos += 12 + chunkLength;
			*cast(uint*)&data[i] = bswap(chunkLength);
			(cast(char[])data[i+4 .. i+8])[] = chunk.type;
			data[i+8 .. i+8+chunk.data.length] = cast(ubyte[])chunk.data;
			*cast(uint*)&data[i+8+chunk.data.length] = bswap(chunk.crc32());
			assert(pos == i+12+chunk.data.length);
		}
		std.file.write(filename, data);
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

struct RGB    { ubyte  r, g, b; }
struct RGB16  { ushort r, g, b; }
struct RGBX   { ubyte  r, g, b, x; }
struct RGBX16 { ushort r, g, b, x; }
struct RGBA   { ubyte  r, g, b, a; }
struct RGBA16 { ushort r, g, b, a; }
struct GA     { ubyte  g, a; }
struct GA16   { ushort g, a; }

private
{
	static assert(RGB.sizeof == 3);
	RGB[2] test;
	static assert(test.sizeof == 6);
}

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
			foreach (x; xmin..xmax+1)
				pixelHR(x*HR, itpl(y1*HR, y2*HR, x, x1, x2), c);
		else
			foreach (y; ymin..ymax+1)
				pixelHR(itpl(x1*HR, x2*HR, y, y1, y2), y*HR, c);
	}

	void fineLine(uint x1, uint y1, uint x2, uint y2, COLOR c, uint step=1)
	{
		auto xmin = min(x1, x2);
		auto xmax = max(x1, x2);
		auto ymin = min(y1, y2);
		auto ymax = max(y1, y2);

		if (xmax-xmin > ymax-ymin)
		{
			uint end = xmax*HR;
			for (uint x=xmin*HR; x<=end; x+=step)
				pixelHR(x, itpl(y1*HR, y2*HR, x, x1*HR, x2*HR), c);
		}
		else
		{
			uint end = ymax*HR;
			for (uint y=ymin*HR; y<=end; y+=step)
				pixelHR(itpl(x1*HR, x2*HR, y, y1*HR, y2*HR), y, c);
		}
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
				string[] fields = structFields!T;

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

template ReplaceType(T, FROM, TO)
{
	static if (is(T == FROM))
		alias TO ReplaceType;
	else
	static if (is(T==struct))
		struct ReplaceType
		{
			static string mixFields()
			{
				string s;
				foreach (field; structFields!T)
					s ~= "ReplaceType!(typeof(" ~ T.stringof ~ ".init." ~ field ~ "), FROM, TO) " ~ field ~ ";\n";
				return s;
			}

			//pragma(msg, mixFields());
			mixin(mixFields());
		}
	else
		static assert(0, "Can't replace " ~ T.stringof);
}

// *****************************************************************************

T itpl(T, U)(T low, T high, U r, U rLow, U rHigh)
{
	return cast(T)(low + (cast(Signed!T)high-cast(Signed!T)low) * (cast(Signed!U)r - cast(Signed!U)rLow) / (cast(Signed!U)rHigh - cast(Signed!U)rLow));
}

private string[] structFields(T)()
{
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
	return fields;
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

	static string mixConvert(T)(string srcVar, string destVar, string convArray)
	{
		static if (is(T==struct))
		{
			string s;
			foreach (field; structFields!T)
				s ~= destVar~"."~field~" = "~convArray~"["~srcVar~"."~field~"];";
			return s;
		}
		else
			return destVar~" = "~convArray~"["~srcVar~"];";
	}

	auto image2pix(COLOR)(in Image!COLOR lumImage)
	{
		alias ReplaceType!(COLOR, LUM_COLOR, PIX_COLOR) COLOR2;
		auto pixImage = Image!COLOR2(lumImage.w, lumImage.h);
		foreach (i, p; lumImage.pixels)
			mixin(mixConvert!COLOR(`p`, `pixImage.pixels[i]`, `lum2pix`));
		return pixImage;
	}
}
