from Compiler.types import personal, sgf2n


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


def _reveal_and_encode_blocks(blocks, party):
    # Keep every block and lane in one SIMD vector through socket write-back.
    flat_blocks = sgf2n.concat(blocks)
    # sgf2n.reveal_to() does not propagate the vector size when creating its
    # mask, so request one independent mask per flattened block and lane.
    secret_mask, clear_mask = sgf2n.get_random_input_mask_for(
        party, size=flat_blocks.size
    )
    revealed = personal(
        party, (flat_blocks + secret_mask).reveal() - clear_mask
    )._v
    return [
        (revealed >> shift) & OUTPUT_VALUE_MASK
        for shift in range(0, GF2N_BITS, OUTPUT_VALUE_BITS)
    ]


def reveal_and_encode_share(share, party):
    """
    Reveal and encode an LRSS share for socket output.

    The share is flattened first. Each GF(2^128) block is emitted in
    block-major, lane-major order as four 32-bit values, least-significant
    part first. Every part fits safely in the signed integer representation
    used by cgf2n.write_to_socket().
    """
    return _reveal_and_encode_blocks(flatten_share(share), party)


def reveal_and_encode_robust_share(share, party):
    """
    Reveal and encode a robust LRSS share for socket output.

    The flat layout is the base LRSS share, two blocks per MAC key, and one
    MAC tag per party.
    """
    return _reveal_and_encode_blocks(flatten_robust_share(share), party)
