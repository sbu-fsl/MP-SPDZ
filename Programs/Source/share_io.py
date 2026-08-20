from Compiler.library import if_
from Compiler.types import cgf2n, personal, sgf2n


GF2N_BITS = 128
OUTPUT_VALUE_BITS = 32
OUTPUT_VALUE_MASK = (1 << OUTPUT_VALUE_BITS) - 1


def read_blocks(player, block_count, size):
    """Read ``block_count`` SIMD field elements supplied by ``player``."""
    return [
        sgf2n.get_input_from(player, size=size)
        for _ in range(block_count)
    ]


def flatten_share(share):
    """Flatten a share made up of nested lists and tuples of field elements."""
    if isinstance(share, (list, tuple)):
        return [
            block
            for component in share
            for block in flatten_share(component)
        ]
    return [share]


def reveal_and_encode_share(share, party):
    """Reveal one share to its owner and encode it for socket output.

    Each GF(2^128) block is emitted as four unsigned 32-bit values, least
    significant part first. Blocks and SIMD lanes retain the order produced by
    ``flatten_share()``.
    """
    blocks = flatten_share(share)
    if not blocks:
        raise ValueError("cannot encode an empty share")

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


def write_shares_to_socket(socket, shares):
    """Reveal and write each share only to the client that owns it."""
    for party, share in enumerate(shares):
        output = reveal_and_encode_share(share, party)

        @if_(party == socket)
        def _():
            cgf2n.write_to_socket(socket, output)
