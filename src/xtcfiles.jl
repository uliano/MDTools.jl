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

mutable struct XtcFile
    file::IOStream
    filename::AbstractString
    mode::AbstractString
    magicints::Vector{Int}
    natoms::Int
    nframes::Int
    steps::Vector{Int32}
    offsets::Vector{UInt64}
    time::Vector{Float32}
    function XtcFile(name::AbstractString, mode::AbstractString)
        file = open(name, mode)
        magicints = [Int(floor(2^i)) for i in 0:1//3:24]
        magicints[1:9] .= 0
        if 'r' in mode || 'a' in mode
            stuff = read_xtc_headers(file)
            XF = new(file, name, mode, magicints, Int(stuff.natoms), length(stuff.offsets), 
            stuff.steps, stuff.offsets, stuff.times)
        else
            XF = new(file, name, mode, magicints, 0, 0, Int32[], UInt64[], Float32[])
        end
        closeme(xf) = close(xf.file)
        finalizer(closeme, XF)
    end
end

function show(io::IO, xf::XtcFile)
    result = "XtcFile(\"$(xf.filename)\" open for \"$(xf.mode)\" "
    result *= "$(xf.nframes) frames with $(xf.natoms) atoms each)"
    println(io, result)
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

function sendbits!(buffer::BitBuffer, bits::Integer, nbits::Integer)
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

function receivebits!(buffer::BitBuffer, nbits::Integer)::Int
    result = 0
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

function sendints!(buffer::BitBuffer, num_of_bits::Integer, sizes::AbstractVector{S}, nums::AbstractVector{T}) where {S,T <: Integer}
    if num_of_bits > 64
        result = Int128(nums[1])
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
        result64 = Int64(nums[1])
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

function receiveints!(buffer::BitBuffer, num_of_bits::Integer, sizes::AbstractVector{S}, nums::AbstractVector{T}) where {S,T <: Integer}
    if num_of_bits > 64 
        result128 = zero(UInt128)
        ones128 = typemax(typeof(result128)) 
        totbits128 = sizeof(result128) * 8
        freebits = totbits128 - num_of_bits
        
        firstbits128 = UInt128(buffer.bits[buffer.index]) << (buffer.offset + 120)
        buffer.index += 1
        num_of_bits -= 8 - buffer.offset
        new_offset = rem(num_of_bits, 8)
        avail_in_last_byte = rem(8 - new_offset, 8)
    
        num_of_bytes = cld(num_of_bits, 8)
        for _ in 1:num_of_bytes 
            result128 <<= 8
            result128 |= buffer.bits[buffer.index]
            buffer.index += 1
        end

        # align left 
        result128 >>>= avail_in_last_byte
    
        # align right
        result128 <<= freebits

        # insert first bits
        result128 |= firstbits128
    
        # do the trick
        shiftbytes = cld(freebits, 8)
        shiftbits = rem(freebits, 8)
        rest128 = ones128 << (8 * shiftbytes)
        part128 = ~ rest128
        result128 = (result128 & rest128) | ((result128 & part128) >>> shiftbits)
        
        result128 = bswap(result128)
    
        nums[3] = rem(result128, sizes[3])
        result64 = div(result128, sizes[3]) % UInt64
        nums[2] = rem(result64, sizes[2])
        nums[1] = div(result64, sizes[2])

    else

        result64 = zero(UInt64)
        ones64 = typemax(typeof(result64)) 
        totbits64 = sizeof(result64) * 8
        freebits = totbits64 - num_of_bits
        
        firstbits64 = UInt64(buffer.bits[buffer.index]) << (buffer.offset + 56)
        buffer.index += 1
        num_of_bits -= 8 - buffer.offset
        new_offset = rem(num_of_bits, 8)
        avail_in_last_byte = rem(8 - new_offset, 8)
    
        num_of_bytes = cld(num_of_bits, 8)
        for _ in 1:num_of_bytes 
            result64 <<= 8
            result64 |= buffer.bits[buffer.index]
            buffer.index += 1
        end

        # align left 
        result64 >>>= avail_in_last_byte
    
        # align right
        result64 <<= freebits

        # insert first bits
        result64 |= firstbits64
    
        # do the trick
        shiftbytes = cld(freebits, 8)
        shiftbits = rem(freebits, 8)
        rest64 = ones64 << (8 * shiftbytes)
        part64 = ~ rest64
        result64 = (result64 & rest64) | ((result64 & part64) >>> shiftbits)
        
        result64 = bswap(result64)
    
        nums[3] = rem(result64, sizes[3])
        result64 = div(result64, sizes[3])
        nums[2] = rem(result64, sizes[2])
        nums[1] = div(result64, sizes[2])

    end
    buffer.offset = new_offset
    if buffer.offset > 0 
        buffer.index -= 1  # we need to reread the last byte
    end
    return nothing
end

function sizeofint(size::Integer)
    nbits = 0
    while size > 0
        size >>= 1
        nbits += 1
    end    
    return nbits
end

function sizeofints(sizes::AbstractVector)
    product::UInt128 = one(Int128)
    nbits = 0
    for i in 1:length(sizes)
        product *= sizes[i]
    end
    while product > 0
        product >>= 1
        nbits += 1
    end
    return nbits
end

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

function read_xtc_box(file::XtcFile, frame::Integer)
    seek(file.file, file.offsets[frame])
    box = Matrix{Float32}(undef, 3, 3)
    for i in 1:3
        for j in 1:3
            read!(file, box[i, j])
        end
    end
    return box
end

function read_xtc_atoms(file::XtcFile, frame::Integer, coords::AbstractMatrix{T}) where {T<: Real}
    local smallnum::Int
    local smaller::Int
    minint = MVector{3, Int32}(undef)
    maxint = MVector{3, Int32}(undef)
    thiscoord = MVector{3, Int}(undef)
    prevcoord = MVector{3, Int}(undef)
    magicints = file.magicints
    FIRSTINDEX = 9
    seek(file.file, file.offsets[frame] + 36) # we don't read box here
    size = ntoh(read(file.file, Int32))
    if size <= 9
        for j in 1:3
            for i in 1:size
                coords[i, j] = ntoh(read(file.file, Float32))
            end
        end
    else
        precision = ntoh(read(file.file, Float32))
        read!(file.file, minint)
        minint .= ntoh.(minint)
        read!(file.file, maxint)
        maxint .= ntoh.(maxint)
        small_idx = ntoh(read(file.file, Int32))
        smallidx = Int(small_idx)
        smaller = magicints[max(FIRSTINDEX, smallidx - 1) + 1] >>> 1 
        smallnum = magicints[smallidx + 1] >>> 1 
        sizesmall = SA[magicints[smallidx + 1], magicints[smallidx + 1], magicints[smallidx + 1]]
        nbytes = ntoh(read(file.file, Int32))
        buffer = BitBuffer(nbytes)
        readbytes!(file.file, buffer.bits, nbytes)
        sizeint = SA[(maxint[1] - minint[1] + 1),
                    (maxint[2] - minint[2] + 1),
                    (maxint[3] - minint[3] + 1)]
        if any( sizeint .> 2^24-1)
            bitsizeint = sizeofint.(sizeint)
            bitsize = 0
        else    
            bitsize = sizeofints(sizeint)
        end
        i = 1 # atom index
        run = 0
        while i <= file.natoms
            if bitsize == 0
                for ii in 1:3
                    thiscoord[ii] = receivebits!(buffer, bitsizeint[ii])
                end
            else
                receiveints!(buffer, bitsize, sizeint, thiscoord)
            end
            for ii in 1:3
                thiscoord[ii] += minint[ii]
            end
            flag = receivebits!(buffer, 1)
            is_smaller = 0
            if flag == 1
                run = receivebits!(buffer, 5)
                is_smaller = run % 3
                run -= is_smaller
                is_smaller -= 1
            end
            if run == 0
                for ii in 1:3
                    coords[ii, i] = thiscoord[ii] / precision 
                end
                i += 1
            else
                for ii in 1:3
                    prevcoord[ii] = thiscoord[ii]
                end
                for k in 1:3:run
                    receiveints!(buffer, smallidx, sizesmall, thiscoord)
                    for ii in 1:3
                        thiscoord[ii] += prevcoord[ii] - smallnum # smallnum ??? WTF!!
                    end
                    if k == 1 # exchange first with second atom WTF!!
                        thiscoord, prevcoord = prevcoord, thiscoord
                        for ii in 1:3
                            coords[ii,i] = prevcoord[ii] / precision
                        end
                        i += 1
                    else
                        for ii in 1:3
                            prevcoord[ii] = thiscoord[ii]
                        end
                    end
                    for ii in 1:3
                        coords[ii,i] = thiscoord[ii] / precision
                    end
                    i += 1
                end
            end
            smallidx += is_smaller;
            if is_smaller < 0
                smallnum = smaller
                if smallidx > FIRSTINDEX 
                    smaller = magicints[smallidx] >>> 1
                else
                    smaller = 0
                end
            elseif is_smaller > 0
                smaller = smallnum
                smallnum = magicints[smallidx + 1] >>> 1
            end
            sizesmall = SA[magicints[smallidx + 1], magicints[smallidx + 1], magicints[smallidx + 1]]
        end
    end
    return nothing      
end

function read_xtc_file(name)
    xtcfile = XtcFile(name, "r")
    coordinates = Array{Float32}(undef, 3, xtcfile.natoms, xtcfile.nframes)
    for frame in 1:xtcfile.nframes
        read_xtc_atoms(xtcfile, frame, view(coordinates, :, :, frame))
    end
    return coordinates
end




#@code_warntype 
#@time read_xtc_file("step5_2.xtc");
# stats=@timed cc=read_xtc_file("md_adam_g0_amber_rep1.xtc");
# println(stats.time)
#cc=read_xtc_file("diala_nowat.xtc");

#read_xtc_file("step5_2.xtc");

# xf = XtcFile("diala_nowat.xtc", "r")
# atoms = Matrix{Float32}(undef, 3, xf.natoms)
# @code_warntype read_xtc_atoms(xf, 1, atoms)
#@time cc=read_xtc_file("md_adam_g0_amber_rep1.xtc");
