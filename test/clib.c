/*____________________________________________________________________________
 |
 | sendbits - encode num into buf using the specified number of bits
 |
 | This routines appends the value of num to the bits already present in
 | the array buf. You need to give it the number of bits to use and you
 | better make sure that this number of bits is enough to hold the value
 | Also num must be positive.
 |
 */

void sendbits(int buf[], int num_of_bits, int num)
{

    unsigned int   cnt, lastbyte;
    int            lastbits;
    unsigned char* cbuf;

    cbuf     = ((unsigned char*)(buf)) + 3 * sizeof(*buf);
    cnt      = (unsigned int)(buf[0]);
    lastbits = buf[1];
    lastbyte = (unsigned int)(buf[2]);
    while (num_of_bits >= 8)
    {
        lastbyte    = (lastbyte << 8) | ((num >> (num_of_bits - 8)) /* & 0xff*/);
        cbuf[cnt++] = lastbyte >> lastbits;
        num_of_bits -= 8;
    }
    if (num_of_bits > 0)
    {
        lastbyte = (lastbyte << num_of_bits) | num;
        lastbits += num_of_bits;
        if (lastbits >= 8)
        {
            lastbits -= 8;
            cbuf[cnt++] = lastbyte >> lastbits;
        }
    }
    buf[0] = cnt;
    buf[1] = lastbits;
    buf[2] = lastbyte;
    if (lastbits > 0)
    {
        cbuf[cnt] = lastbyte << (8 - lastbits);
    }
}

/*_________________________________________________________________________
 |
 | sizeofint - calculate bitsize of an integer
 |
 | return the number of bits needed to store an integer with given max size
 |
 */

int sizeofint(const int size)
{
    int num         = 1;
    int num_of_bits = 0;

    while (size >= num && num_of_bits < 32)
    {
        num_of_bits++;
        num <<= 1;
    }
    return num_of_bits;
}

/*___________________________________________________________________________
 |
 | sizeofints - calculate 'bitsize' of compressed ints
 |
 | given the number of small unsigned integers and the maximum value
 | return the number of bits needed to read or write them with the
 | routines receiveints and sendints. You need this parameter when
 | calling these routines. Note that for many calls I can use
 | the variable 'smallidx' which is exactly the number of bits, and
 | So I don't need to call 'sizeofints for those calls.
 */

int sizeofints(const int num_of_ints, const unsigned int sizes[])
{
    int          i, num;
    int          bytes[32];
    unsigned int num_of_bytes, num_of_bits, bytecnt, tmp;
    num_of_bytes = 1;
    bytes[0]     = 1;
    num_of_bits  = 0;
    for (i = 0; i < num_of_ints; i++)
    {
        tmp = 0;
        for (bytecnt = 0; bytecnt < num_of_bytes; bytecnt++)
        {
            tmp            = bytes[bytecnt] * sizes[i] + tmp;
            bytes[bytecnt] = tmp & 0xff;
            tmp >>= 8;
        }
        while (tmp != 0)
        {
            bytes[bytecnt++] = tmp & 0xff;
            tmp >>= 8;
        }
        num_of_bytes = bytecnt;
    }
    num = 1;
    num_of_bytes--;
    while (bytes[num_of_bytes] >= num)
    {
        num_of_bits++;
        num *= 2;
    }
    return num_of_bits + num_of_bytes * 8;
}

/*____________________________________________________________________________
 |
 | sendints - send a small set of small integers in compressed format
 |
 | this routine is used internally by xdr3dfcoord, to send a set of
 | small integers to the buffer.
 | Multiplication with fixed (specified maximum ) sizes is used to get
 | to one big, multibyte integer. Allthough the routine could be
 | modified to handle sizes bigger than 16777216, or more than just
 | a few integers, this is not done, because the gain in compression
 | isn't worth the effort. Note that overflowing the multiplication
 | or the byte buffer (32 bytes) is unchecked and causes bad results.
 |
 */

void sendints(int          buf[],
                     const int    num_of_ints,
                     const int    num_of_bits,
                     unsigned int sizes[],
                     unsigned int nums[])
{

    int          i, num_of_bytes, bytecnt;
    unsigned int bytes[32], tmp;

    tmp          = nums[0];
    num_of_bytes = 0;
    do
    {
        bytes[num_of_bytes++] = tmp & 0xff;
        tmp >>= 8;
    } while (tmp != 0);

    for (i = 1; i < num_of_ints; i++)
    {
        // if (nums[i] >= sizes[i])
        // {
        //     fprintf(stderr,
        //             "major breakdown in sendints num %u doesn't "
        //             "match size %u\n",
        //             nums[i],
        //             sizes[i]);
        //     exit(1);
        // }
        /* use one step multiply */
        tmp = nums[i];
        for (bytecnt = 0; bytecnt < num_of_bytes; bytecnt++)
        {
            tmp            = bytes[bytecnt] * sizes[i] + tmp;
            bytes[bytecnt] = tmp & 0xff;
            tmp >>= 8;
        }
        while (tmp != 0)
        {
            bytes[bytecnt++] = tmp & 0xff;
            tmp >>= 8;
        }
        num_of_bytes = bytecnt;
    }
    if (num_of_bits >= num_of_bytes * 8)
    {
        for (i = 0; i < num_of_bytes; i++)
        {
            sendbits(buf, 8, bytes[i]);
        }
        sendbits(buf, num_of_bits - num_of_bytes * 8, 0);
    }
    else
    {
        for (i = 0; i < num_of_bytes - 1; i++)
        {
            sendbits(buf, 8, bytes[i]);
        }
        sendbits(buf, num_of_bits - (num_of_bytes - 1) * 8, bytes[i]);
    }
}


/*___________________________________________________________________________
 |
 | receivebits - decode number from buf using specified number of bits
 |
 | extract the number of bits from the array buf and construct an integer
 | from it. Return that value.
 |
 */

int receivebits(int buf[], int num_of_bits)
{

    int            cnt, num, lastbits;
    unsigned int   lastbyte;
    unsigned char* cbuf;
    int            mask = (1 << num_of_bits) - 1;

    cbuf     = (unsigned char*)(buf) + 3 * sizeof(*buf);
    cnt      = buf[0];
    lastbits = (unsigned int)(buf[1]);
    lastbyte = (unsigned int)(buf[2]);

    num = 0;
    while (num_of_bits >= 8)
    {
        lastbyte = (lastbyte << 8) | cbuf[cnt++];
        num |= (lastbyte >> lastbits) << (num_of_bits - 8);
        num_of_bits -= 8;
    }
    if (num_of_bits > 0)
    {
        if (lastbits < num_of_bits)
        {
            lastbits += 8;
            lastbyte = (lastbyte << 8) | cbuf[cnt++];
        }
        lastbits -= num_of_bits;
        num |= (lastbyte >> lastbits) & ((1 << num_of_bits) - 1);
    }
    num &= mask;
    buf[0] = cnt;
    buf[1] = lastbits;
    buf[2] = lastbyte;
    return num;
}

/*____________________________________________________________________________
 |
 | receiveints - decode 'small' integers from the buf array
 |
 | this routine is the inverse from sendints() and decodes the small integers
 | written to buf by calculating the remainder and doing divisions with
 | the given sizes[]. You need to specify the total number of bits to be
 | used from buf in num_of_bits.
 |
 */

void receiveints(int buf[], const int num_of_ints, int num_of_bits, const unsigned int sizes[], int nums[])
{
    int bytes[32];
    int i, j, num_of_bytes, p, num;

    bytes[0] = bytes[1] = bytes[2] = bytes[3] = 0;
    num_of_bytes                              = 0;
    while (num_of_bits > 8)
    {
        bytes[num_of_bytes++] = receivebits(buf, 8);
        num_of_bits -= 8;
    }
    if (num_of_bits > 0)
    {
        bytes[num_of_bytes++] = receivebits(buf, num_of_bits);
    }
    for (i = num_of_ints - 1; i > 0; i--)
    {
        num = 0;
        for (j = num_of_bytes - 1; j >= 0; j--)
        {
            num      = (num << 8) | bytes[j];
            p        = num / sizes[i];
            bytes[j] = p;
            num      = num - p * sizes[i];
        }
        nums[i] = num;
    }
    nums[0] = bytes[0] | (bytes[1] << 8) | (bytes[2] << 16) | (bytes[3] << 24);
}