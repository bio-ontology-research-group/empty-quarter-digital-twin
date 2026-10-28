"""Evaluate the retained KEGG module expression independently of the model's KO table.

Supports the conjunction, alternatives, complexes and optional subunit syntax
used in the pinned M00165 definition. This audits encoding, not enzyme activity.
"""
import re


def split_top(text: str, separators: str) -> list[str]:
    depth, start, parts = 0, 0, []
    for i, char in enumerate(text):
        if char == "(":
            depth += 1
        elif char == ")":
            depth -= 1
        elif depth == 0 and char in separators:
            parts.append(text[start:i])
            start = i+1
        if depth < 0:
            raise ValueError("Unbalanced module expression")
    if depth:
        raise ValueError("Unbalanced module expression")
    parts.append(text[start:])
    return [p.strip() for p in parts if p.strip()]


def satisfied(expression: str, kos: set[str]) -> bool:
    expression = expression.strip()
    for separators, operation in [(" ", all), (",", any), ("+", all)]:
        parts = split_top(expression, separators)
        if len(parts) > 1:
            return operation(satisfied(p, kos) for p in parts)
    if expression.startswith("(") and expression.endswith(")"):
        return satisfied(expression[1:-1], kos)
    parts = split_top(expression, "-")
    if len(parts) > 1:
        return satisfied(parts[0], kos)  # optional subunits do not constrain presence
    if not re.fullmatch(r"K\d{5}", expression):
        raise ValueError(f"Unsupported module expression: {expression}")
    return expression in kos


def module_definition(snapshot: str) -> str:
    return next(line.split(None, 1)[1] for line in snapshot.splitlines()
                if line.startswith("DEFINITION"))
