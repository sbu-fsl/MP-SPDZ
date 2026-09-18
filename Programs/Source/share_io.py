from Compiler.library import get_program, if_
from Compiler.types import Array, cgf2n, personal, sgf2n


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

    # Stage each word in memory to avoid concat's per-lane register expansion.
    # Socket writes interleave vectors lane-first, so four block-major word
    # arrays preserve block -> lane -> word order in a single socket message.
    lane_count = sum(block.size for block in blocks)
    words = [
        Array(lane_count, cgf2n)
        for _ in range(0, GF2N_BITS, OUTPUT_VALUE_BITS)
    ]
    # Large vector loads exceed the optimizer's memory-dependency budget.
    # Protect staging stores/loads, including reuse across calls, explicitly.
    program = get_program()
    program.open_memory_scope(reveal_and_encode_share)
    try:
        offset = 0
        for block in blocks:
            # reveal_to() does not propagate vector size to its mask allocation.
            secret_mask, clear_mask = sgf2n.get_random_input_mask_for(
                party, size=block.size
            )
            revealed = personal(
                party, (block + secret_mask).reveal() - clear_mask
            )._v
            for word, shift in zip(words, range(0, GF2N_BITS, OUTPUT_VALUE_BITS)):
                word.assign_vector((revealed >> shift) & OUTPUT_VALUE_MASK, offset)
            offset += block.size

        output = [word.get_vector() for word in words]
    finally:
        program.close_memory_scope()
    for word in words:
        word.delete()
    return output


def write_shares_to_socket(socket, shares):
    """Reveal and write each share only to the client that owns it."""
    for party, share in enumerate(shares):
        output = reveal_and_encode_share(share, party)

        @if_(party == socket)
        def _():
            cgf2n.write_to_socket(socket, output)
