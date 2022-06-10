include("../src/xtcfiles.jl")
using  Libdl
lib = Libdl.dlopen("./test/libclib.so")
using Test

function buffer2c!(julia::BitBuffer, c::Array{Cint})
    c[1] = 0
    c[2] = 0
    c[3] = 0
    for i in 0:length(julia.bits) ÷ 4
        index = 1 + i * 4
        if index > length(julia.bits) break end
        c[i+4] = Int32(julia.bits[index]) << 0
        index += 1
        if index > length(julia.bits) break end
        c[i+4] |= Int32(julia.bits[index]) << 8
        index += 1
        if index > length(julia.bits) break end
        c[i+4] |= Int32(julia.bits[index]) << 16
        index += 1
        if index > length(julia.bits) break end
        c[i+4] |= Int32(julia.bits[index]) << 24
    end
end

function buffer2j!(c::Array{Cint}, julia::BitBuffer)
    for i in 0:4:((length(c) - 4) * 4)
        int = c[i ÷ 4 + 4]
        len = length(julia.bits)
        julia.bits[i+1] =int & 0xff
        if i+2 > len return nothing end
        int >>>= 8
        julia.bits[i+2] =int & 0xff
        int >>>= 8
        if i+3 > len return nothing end
        julia.bits[i+3] =int & 0xff
        int >>>= 8
        if i+4 > len return nothing end
        julia.bits[i+4] =int & 0xff
    end
    julia.offset = 0
    julia.index = 1
end

function reset!(buffer::BitBuffer)
    buffer.index = 1
    buffer.offset = 0
end

function test_buffer_conversion(n)
    jb1 = BitBuffer(n)
    jb2 = BitBuffer(n)
    len = n ÷ 4 + 3
    if n % 4 > 0 len += 1 end
    cb = Array{Cint}(undef, len)
    for i in 1:length(jb1.bits)
        jb1.bits[i] = rand(1:255) % UInt8
    end
    buffer2c!(jb1, cb)
    buffer2j!(cb, jb2)
    return jb1.bits == jb2.bits[1:length(jb1.bits)]
end

function test_sizeofint(n::Int64)
    sym = Libdl.dlsym(lib, :sizeofint)
    s = sizeofint(n % UInt32)
    r = ccall(sym, Cint, (Cint,), n)
    return r == s
end

function test_sizeofints()
    sym = Libdl.dlsym(lib, :sizeofints)
    sizes = MVector{3, UInt32}(undef)
    result = true
    for i in 2:24
        for j in 2:24
            for k in 2:24
                sizes[1] = UInt32(2^i)
                sizes[2] = UInt32(2^j)
                sizes[3] = UInt32(2^k)
                r = sizeofints(sizes)
                s = ccall(sym, Cuint, (Cint, Ptr{Cuint}), 3, sizes)
                result &= r==s
                if r!=s 
                    println(sizes, " ", r, " ", s)
                end
                sizes[1] = UInt32(2^i - 1)
                sizes[2] = UInt32(2^j)
                sizes[3] = UInt32(2^k)
                r = sizeofints(sizes)
                s = ccall(sym, Cuint, (Cint, Ptr{Cuint}), 3, sizes)
                result &= r==s
                if r!=s 
                    println(sizes, " ", r, " ", s)
                end

                sizes[1] = UInt32(2^i)
                sizes[2] = UInt32(2^j - 1)
                sizes[3] = UInt32(2^k)
                r = sizeofints(sizes)
                s = ccall(sym, Cuint, (Cint, Ptr{Cuint}), 3, sizes)
                result &= r==s
                if r!=s 
                    println(sizes, " ", r, " ", s)
                end

                sizes[1] = UInt32(2^i)
                sizes[2] = UInt32(2^j)
                sizes[3] = UInt32(2^k - 1)
                r = sizeofints(sizes)
                s = ccall(sym, Cuint, (Cint, Ptr{Cuint}), 3, sizes)
                result &= r==s
                if r!=s 
                    println(sizes, " ", r, " ", s)
                end

                sizes[1] = UInt32(2^i)
                sizes[2] = UInt32(2^j - 1)
                sizes[3] = UInt32(2^k - 1)
                r = sizeofints(sizes)
                s = ccall(sym, Cuint, (Cint, Ptr{Cuint}), 3, sizes)
                result &= r==s
                if r!=s 
                    println(sizes, " ", r, " ", s)
                end

                sizes[1] = UInt32(2^i - 1)
                sizes[2] = UInt32(2^j)
                sizes[3] = UInt32(2^k - 1)
                r = sizeofints(sizes)
                s = ccall(sym, Cuint, (Cint, Ptr{Cuint}), 3, sizes)
                result &= r==s
                if r!=s 
                    println(sizes, " ", r, " ", s)
                end

                sizes[1] = UInt32(2^i - 1)
                sizes[2] = UInt32(2^j - 1)
                sizes[3] = UInt32(2^k)
                r = sizeofints(sizes)
                s = ccall(sym, Cuint, (Cint, Ptr{Cuint}), 3, sizes)
                result &= r==s
                if r!=s 
                    println(sizes, " ", r, " ", s)
                end

                sizes[1] = UInt32(2^i - 1)
                sizes[2] = UInt32(2^j - 1)
                sizes[3] = UInt32(2^k - 1)
                r = sizeofints(sizes)
                s = ccall(sym, Cuint, (Cint, Ptr{Cuint}), 3, sizes)
                result &= r==s
                if r!=s 
                    println(sizes, " ", r, " ", s)
                end
            end
        end
    end
    return result
end

function test_send_receive_bits()
    result = true
    rec = Libdl.dlsym(lib, :receivebits)
    send = Libdl.dlsym(lib, :sendbits)
    jbuf = BitBuffer(59)
    cbuf = Array{Cint}(undef, 3+ Int(cld(59, 4)))
    nums = Vector{Int32}(undef, 30)
    jnums = Vector{Int32}(undef, 30)
    cnums = Vector{Int32}(undef, 30)
    sizes = Vector{Int32}(undef, 30)
    for i in 1:30
        nums[i] = Int32(2^i - 1)
        sizes[i] = UInt32(sizeofint(2^i -1))
        sendbits!(jbuf, nums[i], sizes[i])
    end
    buffer2c!(jbuf,cbuf)
    for i in 1:30
        cnums[i] = ccall(rec, Cint, (Ptr{Cint}, Cint), cbuf, sizes[i])
    end
    result &= nums == cnums
    cbuf .= 0
    for i in 1:30
        ccall(send, Cvoid, (Ptr{Cint}, Cint, Cint), cbuf, sizes[i], cnums[i])
    end
    reset!(jbuf)
    for i in 1:30
        jnums[i] = receivebits!(jbuf, sizes[i])
    end
    result &= nums == jnums
    return result
end

function test_send_receive_ints(size)
    result = true
    rec = Libdl.dlsym(lib, :receiveints)
    send = Libdl.dlsym(lib, :sendints)
    sizes = @MVector UInt32[size, size, size]
    isizes = @MVector [size, size, size]
    num_of_bits = sizeofints(sizes)
    start = 0
    stop = size-1
    step = 1
    bytes = cld(size * num_of_bits, 8)
    words = cld(bytes, 4) + 3
    bj=BitBuffer(bytes)
    for i in start:step:stop
        nums = @MVector [i, i, stop - i]
        sendints!(bj, num_of_bits, isizes, nums)
    end
    reset!(bj)
    for i in start:step:stop
        nums =  [i, i, stop - i]
        nums2 = MVector{3,Int}(undef)
        receiveints!(bj, num_of_bits, isizes, nums2)
        result &= nums == nums2
        if nums != nums2 
            println(i, " ", nums, " ", nums2)
        end
    end
    reset!(bj)
    for i in start:step:stop
        nums = @MVector [i, i, stop - i]
        sendints!(bj,num_of_bits, isizes, nums)
    end
    bc = Array{Cint}(undef, words)
    buffer2c!(bj, bc)
    for i in start:step:stop
        nums =  Int32[i, i, stop - i]
        nums2 = Vector{Int32}(undef, 3)
        ccall(rec, Cvoid, (Ptr{Cint}, Cint, Cint, Ptr{Cuint}, Ptr{Cint}), bc, 3, num_of_bits, sizes, nums2)
        result &= nums == nums2
        if nums != nums2 
            println(i, " ", nums, " ", nums2)
        end
    end
    bc .= 0
    for i in start:step:stop
        nums =  UInt32[i, i, stop - i]
        ccall(send, Cvoid, (Ptr{Cint}, Cint, Cint, Ptr{Cuint}, Ptr{Cint}), bc, 3, num_of_bits, sizes, nums)
    end
    buffer2j!(bc, bj)
    reset!(bj)
    for i in start:step:stop
        nums = [i, i, stop - i]
        nums2 = MVector{3,Int}(undef)
        receiveints!(bj, num_of_bits, isizes, nums2)
        result &= nums == nums2
        if nums != nums2 
            println(i, " ", nums, " ", nums2)
        end
    end
    return result
end

function test_send_receive_ints_(nb)
    magicints = [Int(floor(2^i)) for i in 0:1//3:24]
    ll = length(magicints)
    list = []
    for i in 1:ll
        for j in 1:ll
            for k in 1:ll
                bits = sizeofints([magicints[i], magicints[j], magicints[k]])
                if bits == nb
                    push!(list, [magicints[i], magicints[j], magicints[k]])
                end
            end
        end
    end

    result = true
    rec = Libdl.dlsym(lib, :receiveints)
    recb = Libdl.dlsym(lib, :receivebits)
    send = Libdl.dlsym(lib, :sendints)
    sizes = MVector{3,UInt32}(undef)
    bits = 42
    for isizes in list # [[16777215, 16777215, 100]] # [[3329020, 16777215, 511]] # 
        for nbits in 1:8
            for i in 1:3 sizes[i] = convert(UInt32, isizes[i]) end
            num_of_bits = sizeofints(sizes)
            if num_of_bits != 64
                #println("num of bits $num_of_bits")
            end
            bytes = 32
            words = cld(bytes, 4) + 3
            bj=BitBuffer(bytes)
            nums = MVector{3, Int32}(undef)
            for i in 1:3 nums[i] = convert(UInt32, isizes[i] - 1) end
            sendbits!(bj,bits,nbits)
            sendints!(bj, num_of_bits, isizes, nums)
            sendints!(bj, num_of_bits, isizes, nums)
            reset!(bj)
            bc = Array{Cint}(undef, words)
            buffer2c!(bj, bc)
            nums2 = Vector{Int32}(undef, 3)
            bb=ccall(recb, Cint, (Ptr{Cint}, Cint), bc, nbits)
            ccall(rec, Cvoid, (Ptr{Cint}, Cint, Cint, Ptr{Cuint}, Ptr{Cint}), bc, 3, num_of_bits, sizes, nums2)
            ccall(rec, Cvoid, (Ptr{Cint}, Cint, Cint, Ptr{Cuint}, Ptr{Cint}), bc, 3, num_of_bits, sizes, nums2)
            if nums != nums2 
                println(" ", nums, " ", nums2)
                result = false
            end
            nums3 = Vector{Int32}(undef, 3)
            bb=receivebits!(bj,nbits)
            # @code_warntype receiveints!(bj, num_of_bits, isizes, nums3)
            receiveints!(bj, num_of_bits, isizes, nums3)
            receiveints!(bj, num_of_bits, isizes, nums3)
            if nums != nums3 
                println(nbits, "-", nums, "-", nums3)
                result = false
            end
        
            reset!(bj)
        end
    end
    return result
end

@testset "Buffer conversion" begin
    @test test_buffer_conversion(2^24)
    @test test_buffer_conversion(2^24+1)
    @test test_buffer_conversion(2^24+2)
    @test test_buffer_conversion(2^24+3)
end

@testset "sizeofint" begin
    @test test_sizeofint(1)
    @test test_sizeofint(2)
    @test test_sizeofint(3)
    @test test_sizeofint(4)
    @test test_sizeofint(7)
    @test test_sizeofint(8)
    @test test_sizeofint(15)
    @test test_sizeofint(16)
    @test test_sizeofint(31)
    @test test_sizeofint(32)
    @test test_sizeofint(63)
    @test test_sizeofint(64)
    @test test_sizeofint(127)
    @test test_sizeofint(128)
    @test test_sizeofint(255)
    @test test_sizeofint(256)
    @test test_sizeofint(511)
    @test test_sizeofint(512)
    @test test_sizeofint(1023)
    @test test_sizeofint(1024)
    @test test_sizeofint(2047)
    @test test_sizeofint(2048)
    @test test_sizeofint(4095)
    @test test_sizeofint(4096)
    @test test_sizeofint(8191)
    @test test_sizeofint(8192)
    @test test_sizeofint(16383)
    @test test_sizeofint(16384)
    @test test_sizeofint(32767)
    @test test_sizeofint(32768)
    @test test_sizeofint(65535)
    @test test_sizeofint(65536)
    @test test_sizeofint(131071)
    @test test_sizeofint(131072)
    @test test_sizeofint(262143)
    @test test_sizeofint(262144)
    @test test_sizeofint(524287)
    @test test_sizeofint(524288)
    @test test_sizeofint(1048575)
    @test test_sizeofint(1048576)
    @test test_sizeofint(2097151)
    @test test_sizeofint(2097152)
    @test test_sizeofint(4194303)
    @test test_sizeofint(4194304)
    @test test_sizeofint(8388607)
    @test test_sizeofint(8388608)
    @test test_sizeofint(16777215)
    @test test_sizeofint(16777216)
    @test test_sizeofint(33554431)
    @test test_sizeofint(33554432)
    @test test_sizeofint(67108863)
    @test test_sizeofint(67108864)
    @test test_sizeofint(134217727)
    @test test_sizeofint(134217728)
    @test test_sizeofint(268435455)
    @test test_sizeofint(268435456)
    @test test_sizeofint(536870911)
    @test test_sizeofint(536870912)
    @test test_sizeofint(1073741823)
    # in this range the original implementation gives unexpectedly 32, instead of 31,31
    @test_broken test_sizeofint(1073741824)
    @test_broken test_sizeofint(2147483647)
    # original signed implementaion doesn't work for 32 bits as numbers are negative
end

@testset "sizeofints" begin
    @test test_sizeofints()
end

@testset "sendbits! and receivebits!" begin
    @test test_send_receive_bits()
end

@testset "sendints! and receiveints!" for num in 8:24
    # 8 is the lowest limit in magicints
    @test test_send_receive_ints(num)
end

@testset "receiveints!" for nbits in 8:72

    @test test_send_receive_ints_(nbits)
end

Libdl.dlclose(lib);


