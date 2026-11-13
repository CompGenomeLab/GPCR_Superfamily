"""Minimal newick parser: returns root node with .children, .label, .length, .support"""

class Node:
    __slots__ = ("children", "parent", "label", "length", "support", "leaves")
    def __init__(self):
        self.children = []; self.parent = None
        self.label = None; self.length = None; self.support = None
        self.leaves = None

def parse(text):
    text = text.strip()
    i = 0
    root = Node()
    cur = root
    n = len(text)
    while i < n:
        c = text[i]
        if c == "(":
            ch = Node(); ch.parent = cur; cur.children.append(ch); cur = ch; i += 1
        elif c == ",":
            ch = Node(); ch.parent = cur.parent; cur.parent.children.append(ch); cur = ch; i += 1
        elif c == ")":
            cur = cur.parent; i += 1
        elif c == ";":
            break
        else:
            j = i
            while j < n and text[j] not in "(),;":
                j += 1
            tok = text[i:j].strip()
            i = j
            if ":" in tok:
                name, _, ln = tok.rpartition(":")
                try: cur.length = float(ln)
                except ValueError: pass
            else:
                name = tok
            name = name.strip()
            if name:
                if cur.children:
                    cur.support = name          # internal node support value
                else:
                    cur.label = name
    return root

def walk(node):
    yield node
    for c in node.children:
        yield from walk(c)

def leafset(root):
    """populate .leaves (frozenset of leaf labels) bottom-up"""
    order = list(walk(root))
    for nd in reversed(order):
        if not nd.children:
            nd.leaves = frozenset([nd.label])
        else:
            s = set()
            for c in nd.children:
                s |= c.leaves
            nd.leaves = frozenset(s)
    return root
