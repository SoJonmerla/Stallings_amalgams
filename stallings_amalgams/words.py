from collections import deque


def inverse_label(label: str) -> str:
    """
    Return the formal inverse of a generator label.

    Examples
    --------
    >>> inverse_label("a")
    'a^-1'

    >>> inverse_label("a^-1")
    'a'
    """
    if label.endswith("^-1"):
        return label[:-3]

    return label + "^-1"

def inverse_word(word: list[str]) -> list[str]:
    """
    Return the formal inverse of a word.

    The inverse is obtained by reversing the order of the letters and
    replacing every letter by its formal inverse.
    """
    return [
        inverse_label(label) for
    label in reversed(word)
    ]

def free_red(w: list[str]) -> list[str]:
    queue = deque(w.copy())
    red_word = [queue.popleft()]
    while queue:
        letter = queue.popleft()
        if letter == inverse_label(red_word[-1]):
            red_word.pop()
            continue
        red_word.append(letter)
    return red_word

def tree_word(
    vertex: int,
    parent: dict[int, int | None],
    parent_label: dict[int, str]
) -> list[str]:
    """
    Return the word labelling the tree path from the root to vertex.
    """
    word = []
    current = vertex

    while parent[current] is not None:
        word.append(parent_label[current])
        current = parent[current]

    word.reverse()
    return word