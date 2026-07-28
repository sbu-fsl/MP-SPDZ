from Compiler.types import cint, regint, sgf2n


GF2N_BITS = 128
OUTPUT_VALUE_BITS = 32
OUTPUT_VALUE_MASK = (1 << OUTPUT_VALUE_BITS) - 1


def read_share(player, num_parties, source_length, size):
    """Read one SIMD batch of LRSS shares in source|ct|seed|masks order."""
    source = [
        sgf2n.get_input_from(player, size=size)
        for _ in range(source_length)
    ]
    ciphertext = sgf2n.get_input_from(player, size=size)
    seed_shares = [
        sgf2n.get_input_from(player, size=size)
        for _ in range(source_length)
    ]
    mask_shares = [
        sgf2n.get_input_from(player, size=size)
        for _ in range(num_parties)
    ]
    return source, ciphertext, seed_shares, mask_shares


def flatten_share(share):
    """Flatten an LRSS share in the canonical source|ct|seed|masks order."""
    source, ciphertext, seed_shares, mask_shares = share
    return [*source, ciphertext, *seed_shares, *mask_shares]


def reveal_and_encode_share(share, party, size):
    """
    Reveal and encode an LRSS share for socket output.

    The share is flattened first. Each GF(2^128) block is emitted in
    block-major, lane-major order as four 32-bit values, least-significant
    part first. This keeps socket values safely inside the arithmetic
    clear-value domain.
    """
    output = []
    for block in flatten_share(share):
        revealed = block.reveal_to(party)
        for lane in range(size):
            clear_value = revealed[lane]._v
            output.extend(
                cint(regint((clear_value >> shift) & OUTPUT_VALUE_MASK))
                for shift in range(0, GF2N_BITS, OUTPUT_VALUE_BITS)
            )
    return output
