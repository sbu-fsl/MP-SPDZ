from Compiler.types import cint, regint, sgf2n


GF2N_BITS = 128
OUTPUT_VALUE_BITS = 32
OUTPUT_VALUE_MASK = (1 << OUTPUT_VALUE_BITS) - 1


def read_share(player, num_parties, source_length, size):
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
    source, ciphertext, seed_shares, mask_shares = share
    return [*source, ciphertext, *seed_shares, *mask_shares]


def read_robust_share(player, num_parties, source_length, size):
    lrss_share = read_share(player, num_parties, source_length, size)
    mac_keys = [
        (
            sgf2n.get_input_from(player, size=size),
            sgf2n.get_input_from(player, size=size),
        )
        for _ in range(num_parties)
    ]
    mac_tags = [
        sgf2n.get_input_from(player, size=size)
        for _ in range(num_parties)
    ]
    return lrss_share, mac_keys, mac_tags


def flatten_robust_share(share):
    lrss_share, mac_keys, mac_tags = share
    mac_key_blocks = []
    for key in mac_keys:
        mac_key_blocks.extend(key)
    return [*flatten_share(lrss_share), *mac_key_blocks, *mac_tags]


def _reveal_and_encode_blocks(blocks, party, size):
    output = []
    for block in blocks:
        revealed = block.reveal_to(party)
        for lane in range(size):
            clear_value = revealed[lane]._v
            output.extend(
                cint(regint((clear_value >> shift) & OUTPUT_VALUE_MASK))
                for shift in range(0, GF2N_BITS, OUTPUT_VALUE_BITS)
            )
    return output


def reveal_and_encode_share(share, party, size):
    """
    Reveal and encode an LRSS share for socket output.

    The share is flattened first. Each GF(2^128) block is emitted in
    block-major, lane-major order as four 32-bit values, least-significant
    part first. This keeps socket values safely inside the arithmetic
    clear-value domain.
    """
    return _reveal_and_encode_blocks(flatten_share(share), party, size)


def reveal_and_encode_robust_share(share, party, size):
    """
    Reveal and encode a robust LRSS share for socket output.

    The flat layout is the base LRSS share, two blocks per MAC key, and one
    MAC tag per party.
    """
    return _reveal_and_encode_blocks(flatten_robust_share(share), party, size)
