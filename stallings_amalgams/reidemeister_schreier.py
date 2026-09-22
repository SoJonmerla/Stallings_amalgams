"""
Presentation simplification following

    G. Havas,
    "A Reidemeister-Schreier program",
    Lecture Notes in Mathematics 372 (1974), §2.4 - §2.6.

The top-level entry point is `simplify_presentation`.
"""

from __future__ import annotations

from .presentation import Presentation
from .words import inverse_label, inverse_word, free_reduce


# ---------------------------------------------------------------------------
# Letter / word helpers
# ---------------------------------------------------------------------------

def split_letter(letter: str) -> tuple[str, int]:
    """Return (generator, sign) with sign = +1 for ``g``, -1 for ``g^-1``."""
    if letter.endswith("^-1"):
        return letter[:-3], -1
    return letter, 1



def cyclically_reduce(word: list[str]) -> list[str]:
    """Free reduce, then cancel inverse pairs wrapping around the ends."""
    word = free_reduce(word)
    while word and word[0] == inverse_label(word[-1]):
        word = word[1:-1]
    return word


def cyclic_rotations(word: list[str]):
    n = len(word)
    for i in range(n):
        yield word[i:] + word[:i]


def _word_sort_key(word: list[str], gen_order: dict[str, int]) -> tuple:
    """
    Key for the ordering of Havas §2.4.  A letter ``g^e`` is smaller than a
    letter ``h^d`` when ``g < h``, or when ``g == h`` and ``e = +1``,
    ``d = -1``.  Two words are compared lexicographically letter by letter,
    with shorter prefixes coming first.
    """
    return tuple(
        (gen_order[g], 0 if s == 1 else 1)
        for g, s in (split_letter(l) for l in word)
    )


def canonical_relator(word: list[str], gen_order: dict[str, int]) -> list[str]:
    """
    Canonical form of the conjugacy class of ``word`` in the free group.

    The canonical representative is the least of all cyclic rotations of
    ``word`` and of its formal inverse, with respect to the ordering
    defined in Havas §2.4.
    """
    word = cyclically_reduce(word)
    if not word:
        return []
    inv = inverse_word(word)
    candidates = list(cyclic_rotations(word)) + list(cyclic_rotations(inv))
    return list(min(candidates, key=lambda w: _word_sort_key(w, gen_order)))


# ---------------------------------------------------------------------------
# Substitution
# ---------------------------------------------------------------------------

def _substitute(word: list[str], values: dict[str, list[str]]) -> list[str]:
    """Replace each generator in ``values`` by its value, everywhere in ``word``."""
    out: list[str] = []
    for letter in word:
        g, s = split_letter(letter)
        if g in values:
            v = values[g]
            if s == -1:
                v = inverse_word(v)
            out.extend(v)
        else:
            out.append(letter)
    return out


# ---------------------------------------------------------------------------
# §2.5 - Eliminate redundant generators
# ---------------------------------------------------------------------------

def _eliminate_redundant_generators(
    generators: list[str],
    relators: list[list[str]],
) -> tuple[list[str], list[list[str]], dict[str, list[str]]]:
    """
    Repeatedly eliminate generators that occur exactly once in a relator.

    This follows the *spirit* of Havas §2.5.1 (technique 1): relators are
    examined in order of increasing length, so that the eliminated generator
    tends to acquire a short equivalent.  When a redundant generator is
    found, its defining relator is removed, the generator is marked
    redundant, and its value is substituted into every remaining relator.
    The process repeats until no redundancy remains.

    One generator is eliminated per inner
    pass; in the terminology of Havas this is closer to technique 2, but
    the outcome on small inputs is the same as technique 1.

    Returns
    -------
    generators : list[str]
        Surviving generators.
    relators : list[list[str]]
        Relators with all redundancies substituted out.
    eliminations : dict[str, list[str]]
        Eliminated generator -> its value as a word in surviving generators.
    """
    generators = list(generators)
    relators = [list(r) for r in relators]
    eliminations: dict[str, list[str]] = {}

    while True:
        # Shorter relators first.
        order = sorted(range(len(relators)), key=lambda i: len(relators[i]))
        progress = False

        for i in order:
            r = relators[i]

            # Count surviving generators appearing in this relator.
            counts: dict[str, int] = {}
            for letter in r:
                g, _ = split_letter(letter)
                if g in eliminations:
                    continue  # already substituted out, should not happen
                counts[g] = counts.get(g, 0) + 1

            redundant: str | None = None
            for g, c in counts.items():
                if c == 1:
                    redundant = g
                    break
            if redundant is None:
                continue

            # Locate the redundant letter in the relator.
            for j, letter in enumerate(r):
                g, s = split_letter(letter)
                if g != redundant:
                    continue
                rest = r[:j] + r[j + 1:]
                value = inverse_word(rest) if s == 1 else rest
                value = free_reduce(value)

                eliminations[redundant] = value
                generators.remove(redundant)
                relators.pop(i)

                # Substitute into every remaining relator.
                new_relators: list[list[str]] = []
                for r2 in relators:
                    r2 = _substitute(r2, {redundant: value})
                    r2 = cyclically_reduce(r2)
                    if r2:
                        new_relators.append(r2)
                relators = new_relators

                progress = True
                break

            if progress:
                break

        if not progress:
            break

    return generators, relators, eliminations


# ---------------------------------------------------------------------------
# §2.6 - Simplify using known orders
# ---------------------------------------------------------------------------

def _discover_orders(relators: list[list[str]]) -> dict[str, int]:
    """Return ``{g: n}`` for every relator of the form ``g^n`` with ``n >= 2``."""
    orders: dict[str, int] = {}
    for r in relators:
        if not r:
            continue
        gens: set[str] = set()
        total = 0
        for letter in r:
            g, s = split_letter(letter)
            gens.add(g)
            total += s
        if len(gens) == 1:
            n = abs(total)
            if n >= 2:
                g = next(iter(gens))
                if g not in orders or n < orders[g]:
                    orders[g] = n
    return orders


def _reduce_powers(word: list[str], orders: dict[str, int]) -> list[str]:
    """Replace runs of a generator of known order by an equivalent shorter run."""
    out: list[str] = []
    i = 0
    n = len(word)
    while i < n:
        g, _ = split_letter(word[i])
        if g not in orders:
            out.append(word[i])
            i += 1
            continue
        order = orders[g]
        exponent = 0
        j = i
        while j < n:
            gj, sj = split_letter(word[j])
            if gj != g:
                break
            exponent += sj
            j += 1
        e = exponent % order
        # Prefer the representative of smallest absolute value.
        if e > order // 2:
            e -= order
        elif e < -(order // 2):
            e += order
        if e > 0:
            out.extend([g] * e)
        elif e < 0:
            out.extend([f"{g}^-1"] * (-e))
        i = j
    return out


def _simplify_by_orders(
    relators: list[list[str]],
) -> tuple[list[list[str]], bool]:
    """
    Havas §2.6.1: use known orders to shorten relators.

    Returns ``(new_relators, changed)``.
    """
    orders: dict[str, int] = {}
    defining: set[int] = set()
    for i, r in enumerate(relators):
        if not r:
            continue
        gens: set[str] = set()
        total = 0
        for letter in r:
            g, s = split_letter(letter)
            gens.add(g)
            total += s
        if len(gens) == 1:
            n = abs(total)
            if n >= 2:
                g = next(iter(gens))
                if g not in orders:
                    orders[g] = n
                    defining.add(i)

    if not orders:
        return relators, False

    changed = False
    new_relators = []
    for i, r in enumerate(relators):
        if i in defining:
            new_relators.append(r)
            continue
        new_r = _reduce_powers(r, orders)
        new_r = cyclically_reduce(new_r)
        if new_r != r:
            changed = True
        if new_r:
            new_relators.append(new_r)
    return new_relators, changed

# ---------------------------------------------------------------------------
# Top-level
# ---------------------------------------------------------------------------

def simplify_presentation(
    P: Presentation,
    max_iterations: int = 1000,
) -> Presentation:
    """
    Simplify a finite presentation following Havas §2.4 - §2.6.

    Steps
    -----
    1. Canonicalise every relator and drop duplicates.
    2. Repeatedly apply §2.5 (eliminate redundant generators) and
       §2.6 (shorten relators using known orders) until neither changes
       the presentation.

    Parameters
    ----------
    P : Presentation
        The presentation to simplify.
    max_iterations : int
        Safety bound on the outer loop.

    Returns
    -------
    Presentation
        The simplified presentation, with relators in canonical form
        and only the surviving generators listed.
    """
    generators = list(P.generators)
    gen_order = {g: i for i, g in enumerate(generators)}

    # Canonicalise and deduplicate the input relators.
    seen: dict[tuple, list[str]] = {}
    for r in P.relators:
        c = canonical_relator(r, gen_order)
        if c:
            seen[tuple(c)] = c
    relators = list(seen.values())

    for _ in range(max_iterations):
        generators, relators, eliminations = _eliminate_redundant_generators(
            generators, relators
        )

        relators, changed_26 = _simplify_by_orders(relators)

        # Canonicalise and deduplicate again.
        gen_order = {g: i for i, g in enumerate(generators)}
        seen = {}
        for r in relators:
            c = canonical_relator(r, gen_order)
            if c:
                seen[tuple(c)] = c
        new_relators = list(seen.values())

        if not eliminations and not changed_26:
            relators = new_relators
            break
        relators = new_relators

    # Drop entries for eliminated generators, leaving the rest untouched.
    # (Generator words are words in the *ambient* group, so they are not
    # affected by internal presentation eliminations.)
    new_generator_words = {
        g: w for g, w in P.generator_words.items() if g in generators
    }

    return Presentation(generators, relators, new_generator_words)