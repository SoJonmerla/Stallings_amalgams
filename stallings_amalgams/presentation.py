from dataclasses import dataclass, field


@dataclass
class Presentation:
    """
    A finite group presentation.
    """

    generators: list[str]
    relators: list[list[str]]
    generator_words: dict[str, list[str]] = field(
        default_factory=dict
    )

    def __str__(self) -> str:
        generators = ", ".join(self.generators)

        if self.relators:
            relators = ", ".join(
                self.format_word(relator)
                for relator in self.relators
            )
        else:
            relators = ""

        return f"< {generators} | {relators} >"

    @staticmethod
    def format_word(word: list[str]) -> str:
        """Convert a word represented by labels into readable text."""
        if not word:
            return "1"

        return " ".join(word)