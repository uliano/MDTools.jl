import Base: show
using StaticArrays

mutable struct BitBuffer
    bits::Vector{UInt8}  # the bytes
    index::Int      # last byte 
    offset::Int     # last bit

    function BitBuffer(data::Vector{UInt8})
        new(data, 1, 0)
    end 
    function BitBuffer(n)
        new(Vector{UInt8}(undef, n), 1, 0)
    end
end

function show(io::IO, buffer::BitBuffer) 
    result = ""
    if buffer.index == 0 return nothing end
    nbits = length(buffer.bits) * 8 - (8 - buffer.offset)
    if buffer.offset == 0 nbits += 8 end
    index = 1
    while nbits >= 8
        result *= bitstring(buffer.bits[index])
        index += 1
        nbits -=8
    end
    if nbits > 0
        result *= bitstring(buffer.bits[index])[1:nbits]
    end
    println(io, result, " - ", length(result), " bit(s).")
end

function send!(buffer::BitBuffer, bits::Integer, nbits::Integer)
    if buffer.offset > 0
        avail_in_first_byte = 8 - buffer.offset
        bits_in_first_byte = min(nbits, avail_in_first_byte)
        byte = ((bits >> (nbits - bits_in_first_byte) % UInt8) & (0xFF >> buffer.offset))
        shift = avail_in_first_byte - bits_in_first_byte
        buffer.bits[buffer.index] |= byte << shift
        buffer.offset = mod(buffer.offset + bits_in_first_byte, 8)
        if buffer.offset == 0 buffer.index += 1 end
        nbits -= bits_in_first_byte
        if nbits <= 0 return nothing end
    end
    while nbits >= 8
        byte = (bits >> (nbits - 8)) % UInt8
        buffer.bits[buffer.index] = byte
        buffer.index += 1
        nbits -= 8
    end
    if nbits > 0
        mask = 0xFF >> (8 - nbits)
        byte = ((bits & mask) % UInt8) << (8 - nbits)
        buffer.bits[buffer.index] = byte
        buffer.offset = nbits
    end
    return nothing
end

function receive!(buffer::BitBuffer, nbits::Integer, returntype::DataType=Int32) 
    result = zero(returntype)
    if buffer.offset > 0
        avail_in_first_byte = 8 - buffer.offset
        bits_in_first_byte = min(nbits, avail_in_first_byte)
        shift = avail_in_first_byte - bits_in_first_byte
        result |= (buffer.bits[buffer.index] & (0xFF >> buffer.offset)) >> shift
        buffer.offset = mod(buffer.offset + bits_in_first_byte, 8)
        if buffer.offset == 0 buffer.index += 1 end
        nbits -= bits_in_first_byte
        if nbits <= 0 return result end
    end
    while nbits >= 8
        result <<= 8
        result |= buffer.bits[buffer.index]
        buffer.index += 1
        nbits -= 8
    end
    if nbits > 0
        result <<= nbits
        mask = 0xFF << (8 - nbits)
        byte = (buffer.bits[buffer.index] & mask) >> (8 - nbits)
        result |= byte
        buffer.offset = nbits
    end
    return result
end

function sendints!(buffer::BitBuffer, num_of_bits::Integer, sizes::MVector{3, UInt32}, nums::MVector{3, UInt32}) 
    if num_of_bits > 64
        result = UInt128(nums[1])
        result *= sizes[2]
        result += nums[2]
        result *= sizes[3]
        result += nums[3]

        avail_in_first_byte = mod(8 - buffer.offset, 8)
        freebits = sizeof(result) * 8 - num_of_bits

        # send first bits
        if buffer.offset > 0  
            firstbits = ((result & (0xFF << buffer.offset)) % UInt8) >>> buffer.offset
            buffer.bits[buffer.index] |= firstbits
            buffer.index += 1
            num_of_bits -= avail_in_first_byte
        end

        # first swap needed to smoothly shift between bytes
        result = bswap(result)

        # this correction is needed to get rid of leading zeoroes 
        # that the bitswap sneaked in between
        shiftbytes = cld(freebits, 8)
        shiftbits = rem(freebits, 8)   
        rest = typemax(typeof(result)) << (8 * shiftbytes)
        part = ~ rest
        value = (result & part) << shiftbits
        result = (result & rest) | value

        # this is the motivation for the first bitswap
        result <<= avail_in_first_byte

        # second swap needed to facilitate sending bytes in little endian 
        result = bswap(result)

        num_of_bytes = cld(num_of_bits, 8)
        for _ in 1:num_of_bytes
            byte = result % UInt8
            buffer.bits[buffer.index] = byte
            buffer.index += 1
            result >>>= 8
        end
    else  # num_of_bits <= 64
        result64 = UInt64(nums[1])
        result64 *= sizes[2]
        result64 += nums[2]
        result64 *= sizes[3]
        result64 += nums[3]

        avail_in_first_byte = mod(8 - buffer.offset, 8)
        freebits = sizeof(result64) * 8 - num_of_bits

        # send first bits
        if buffer.offset > 0  
            firstbits = ((result64 & (0xFF << buffer.offset)) % UInt8) >>> buffer.offset
            buffer.bits[buffer.index] |= firstbits
            buffer.index += 1
            num_of_bits -= avail_in_first_byte
        end

        # first swap needed to smoothly shift between bytes
        result64 = bswap(result64)

        # this correction is needed to get rid of leading zeoroes 
        # that the bitswap sneaked in between
        shiftbytes = cld(freebits, 8)
        shiftbits = rem(freebits, 8)   
        rest64 = typemax(typeof(result64)) << (8 * shiftbytes)
        part64 = ~ rest64
        value64 = (result64 & part64) << shiftbits
        result64 = (result64 & rest64) | value64

        # this is the motivation for the first bitswap
        result64 <<= avail_in_first_byte

        # second swap needed to facilitate sending bytes in little endian 
        result64 = bswap(result64)

        num_of_bytes = cld(num_of_bits, 8)
        for _ in 1:num_of_bytes
            byte = result64 % UInt8
            buffer.bits[buffer.index] = byte
            buffer.index += 1
            result64 >>>= 8
        end
    end
    buffer.offset = mod(num_of_bits, 8)
    if buffer.offset > 0 
        buffer.index -= 1 # last byte is not complete
    end
    return nothing
end

function receiveints!(buffer::BitBuffer, num_of_bits::Integer, sizes::MVector{3, UInt32}, nums::MVector{3, Int32})
    if num_of_bits > 64
        result = zero(UInt128)

    ones = typemax(typeof(result))
    totbits = sizeof(result) * 8
    freebits = totbits - num_of_bits
    avail_in_first_byte = rem(8 - buffer.offset, 8)

    # receive first bits 
    if buffer.offset > 0   
        result |= buffer.bits[buffer.index] & (0xff >>> buffer.offset)
        buffer.index += 1
        num_of_bits -= avail_in_first_byte
    end

    num_of_bytes = cld(num_of_bits, 8)
    for _ in 1:num_of_bytes 
        result <<= 8
        result |= buffer.bits[buffer.index]
        buffer.index += 1
    end

    # align right
    result <<= ((sizeof(result) - num_of_bytes) * 8 - avail_in_first_byte) 

    # clear trailing bits 
    result &= (ones << freebits)

    # do the trick
    shiftbytes = cld(freebits, 8)
    shiftbits = rem(freebits, 8)
    rest = ones << (8 * shiftbytes)
    part = ~ rest
    result = (result & rest) | ((result & part) >>> shiftbits)
    
    result = bswap(result)

    nums[3] = rem(result, sizes[3]) % Int32
    result = div(result, sizes[3]) % UInt64
    nums[2] = rem(result, sizes[2]) % Int32
    nums[1] = div(result, sizes[2]) % Int32

    else
        result64 = zero(UInt64)
        ones64 = typemax(typeof(result64))
        totbits64 = sizeof(result64) * 8
        freebits = totbits64 - num_of_bits
        avail_in_first_byte = rem(8 - buffer.offset, 8)
    
        # receive first bits 
        if buffer.offset > 0   
            result64 |= buffer.bits[buffer.index] & (0xff >>> buffer.offset)
            buffer.index += 1
            num_of_bits -= avail_in_first_byte
        end
    
        num_of_bytes = cld(num_of_bits, 8)
        for _ in 1:num_of_bytes 
            result64 <<= 8
            result64 |= buffer.bits[buffer.index]
            buffer.index += 1
        end
    
        # align right
        result64 <<= ((sizeof(result64) - num_of_bytes) * 8 - avail_in_first_byte) 
    
        # clear trailing bits 
        result64 &= (ones64 << freebits)
    
        # do the trick
        shiftbytes = cld(freebits, 8)
        shiftbits = rem(freebits, 8)
        rest64 = ones64 << (8 * shiftbytes)
        part64 = ~ rest64
        result64 = (result64 & rest64) | ((result64 & part64) >>> shiftbits)
        
        result64 = bswap(result64)
    
        nums[3] = rem(result64, sizes[3]) % Int32
        result64 = div(result64, sizes[3])
        nums[2] = rem(result64, sizes[2]) % Int32
        nums[1] = div(result64, sizes[2]) % Int32

    end

    buffer.offset = rem(num_of_bits, 8)
    if buffer.offset > 0 
        buffer.index -= 1  # we need to reread the last byte
    end
    return nothing
end

function sizeofint(num::Unsigned)
    nbits:: UInt32 = 0
    while num > 0
        num >>= 1
        nbits += 1
    end
    return nbits
end

function sizeofints(sizes::MVector{3, UInt32})
    product::UInt128 = one(UInt128)
    nbits::UInt32 = 0
    for i in 1:length(sizes)
        product *= UInt128(sizes[i])
    end
    while product > 0
        product >>= 1
        nbits += 1
    end
    return nbits
end



const XDR_INT_SIZE = 4
const maxAbsoluteInt = prevfloat(convert(Float32,typemax(Int32)))

magicints = [UInt32(floor(2^i)) for i in 1//3:1//3:24]
magicints[1:8] .= 0

# magicints = Int32[
#     0,        0,        0,       0,       0,       0,       0,       0,       0,       8,
#     10,       12,       16,      20,      25,      32,      40,      50,      64,      80,
#     101,      128,      161,     203,     256,     322,     406,     512,     645,     812,
#     1024,     1290,     1625,    2048,    2580,    3250,    4096,    5060,    6501,    8192,
#     10321,    13003,    16384,   20642,   26007,   32768,   41285,   52015,   65536,   82570,
#     104031,   131072,   165140,  208063,  262144,  330280,  416127,  524287,  660561,  832255,
#     1048576,  1321122,  1664510, 2097152, 2642245, 3329021, 4194304, 5284491, 6658042, 8388607,
#     10568983, 13316085, 16777216
# ]

const FIRSTINDEX = 9
const LASTINDEX = length(magicints)

function read_xtc_headers(file)
    times = Float32[]
    steps = Int32[]
    offsets = UInt64[]
    local natoms
    while ! eof(file)
        natoms, step, time, offset = read_xtc_header(file, skip_to_next=true)
        push!(times, time)
        push!(steps, step)
        push!(offsets, offset)
    end
    return (; natoms, steps, times, offsets)
end

function read_xtc_header(file; skip_to_next=false)
    magic = ntoh(read(file, Int32))
    if magic != 1995 error("Wrong magic number in xtc file.") end
    natoms = ntoh(read(file, Int32))
    step = ntoh(read(file, Int32))
    time = ntoh(read(file, Float32))
    offset = position(file)

    if natoms <=9
        nskip = 4 * (10 + 3 * natoms)
    else
        skip(file, 18 * 4)
        nbytes = ntoh(read(file, Int32))
        nskip = div(nbytes, 4) * 4
        if nbytes % 4 > 0 nskip+= 4 end
    end
    if skip_to_next 
        skip(file, nskip) 
    else
        seek(file, offset)
    end
    return (; natoms, step, time, offset)
end

function read_xtc_data(file, offset)
    
    box = Matrix{Float32}(undef, 3, 3)
    minint = Vector{Int32}(undef, 3)
    maxint = Vector{Int32}(undef, 3)

    seek(file, offset)
    for i in 1:3
        for j in 1:3
            box[i, j] = ntoh(read(file, Float32))
        end
    end
    natoms = ntoh(read(file, Int32))
    coords = Matrix{Float32}(undef, natoms, 3)
    if natoms <= 9
        for j in 1:3
            for i in 1:natoms
                fcoords[i, j] = ntoh(read(file, Float32))
            end
        end
    else
        precision = ntoh(read(file, Float32)) 
        if precision == 0
            # remember to throw error
        end
        read!(file, minint)
        minint .= ntoh.(minint)
        read!(file, maxint)
        maxint .= ntoh.(maxint)
        smallidx = ntoh(read(file, Int32))
        nbytes = ntoh(read(file, Int32))
        bytes = Array{UInt8}(undef, nbytes)
        read!(file, bytes)
        bits = BitBuffer(bytes)

    end
    return (; coords, box)
end


function read_xtc_frame(file)
    box = MMatrix{3, 3, Float32}(undef)
    minint = MVector{3, Int32}(undef)
    maxint = MVector{3, Int32}(undef)
    magic = ntoh(read(file, Int32))
    natoms = ntoh(read(file, Int32))
    step = ntoh(read(file, Int32))
    time = ntoh(read(file, Float32))
    read!(file, box)
    box .= ntoh.(box)
    box = box'
    size = ntoh(read(file, Int32))
    precision = ntoh(read(file, Float32))
    read!(file, minint)
    minint .= ntoh.(minint)
    read!(file, maxint)
    maxint .= ntoh.(maxint)
    smallidx = ntoh(read(file, Int32))
    smaller = magicints[max(FIRSTINDEX, smallidx - 1)] ÷ 2 
    smallnum = magicints[smallidx] ÷ 2 
    sizesmall .= magicints[smallidx]
    println(typeof(smaller))
    nbytes = ntoh(read(file, Int32))
    nskip = cld(nbytes, 4) * 4
    println("smallindex $smallidx")
    #coords = Vector{MVector{3, Float32}}(undef, natoms)
    buffer = BitBuffer(nbytes)
    read!(file, buffer.bits)
    skip(file, nskip-nbytes) 

    sizeint = (maxint .- minint) .% UInt32
    if any( sizeint .> 0xffffff)
        bitsizeint = sizeofint.(sizeint)
        bitsize = 0
    else    
        bitsize = sizeofints(sizeint)
    end
    coords = Matrix{Int32}(undef, 3, natoms)
    prevcoord = MVector{3, Int32}(undef)
    i = 1
    while i <= natoms
        if bitsize == 0
            # ensure order, you never know what the compiler does
            coords[1, i] = receive!(buffer, bitsizeint[1])
            coords[2, i] = receive!(buffer, bitsizeint[2])
            coords[3, i] = receive!(buffer, bitsizeint[3])
        else
            receiveints!(buffer, bitsize, sizeint, ucoord[:])
        end
        coords[:, i] += minint
        prevcoord = coords[:, i]
        i += 1

        flag = receive!(buffer, 1)
        is_smaller = 0
        if flag == 1
            run = receive!(buffer, 5)
            is_smaller = run % 3
            run -= is_smaller
            is_smaller -= 1
        end
        if run > 0
 

        end




    end

    println(bitsize)


    return box, coords      
end

#println(readdir())


xtcfile = open("step5_2.xtc", "r")

fr = read_xtc_headers(xtcfile)
seek(xtcfile, 0)
read_xtc_frame(xtcfile)