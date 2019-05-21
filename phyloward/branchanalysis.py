#!/usr/bin/env python3

# Modified by JaeHeung Han <hylowaker@gmail.com>
# Original code by Yeong Ouk Kim <mrcorndog88@gmail.com>
# Copyright (c) 2017-2018 Yeong Ouk Kim and JaeHeung Han.
import sys

try:
    from typing import Any, Union, List, Dict, Optional, IO
except ImportError:
    print("# ERROR: typing module not found. \n If you are using old version of python, see here 'https://pypi.org/project/typing/'", file=sys.stderr)
    pass


class Node:
    def __init__(self, data: Any=None, name: str=''):
        self.data = data
        self.name = name
        self.distance = 0.
        self.distance_from_root = 0.
        self.parent = None  # type: Optional[Node]
        self.children = []  # type: List[Node]
        self._number_of_descendants = 0
        self._height = 0
        self.level = 0
        self.length = 0.

    @property
    def height(self):
        return self._height

    @property
    def number_of_descendants(self):
        return self._number_of_descendants

    def is_root(self):
        return self.parent is None

    def is_leaf(self):
        return len(self.children) == 0

    def add_child(self, child: 'Node', *, index: int=None, distance: float=None):
        if not child.is_root():
            return False

        if index is None:
            index = len(self.children)
        if distance is not None:
            child.distance = distance
        new_descendants = 1 + child.number_of_descendants
        current, last = self, child

        while True:
            current._number_of_descendants += new_descendants

            new_height = 1 + last._height
            if new_height > current._height:
                current._height = new_height

            new_length = last.length + last.distance
            if new_length > current.length:
                current.length = new_length

            current, last = current.parent, current
            if current is None:
                break

        child.update_level(self.level + 1)
        child.update_distance_from_root(self.distance_from_root + child.distance)
        child.parent = self
        self.children.insert(index, child)
        return True

    def remove_child(self, index: int):
        if not 0 < index < len(self.children):
            return None

        child = self.children[index]
        old_descendants = 1 + child.number_of_descendants
        current = self

        while True:
            current._number_of_descendants -= old_descendants

            max_height = 0
            max_length = 0.
            for i, child in enumerate(self.children):
                if i == index:
                    continue
                if max_height < child.height:
                    max_height = child.height
                if max_length < child.length + child.distance:
                    max_length = child.length + child.distance
            current._height = max_height + 1
            current.length = max_length

            current = current.parent
            if current is None:
                break

    def update_level(self, level: int):
        self.level += level
        for child in self.children:
            child.update_level(level)

    def update_distance_from_root(self, d: float):
        self.distance_from_root += d
        for child in self.children:
            child.update_distance_from_root(d)

    @classmethod
    def calculate_length(cls, current: 'Node'):
        if current.is_leaf():
            return 0.

        max_length = -1 / 0.
        for child in current.children:
            length = child.distance + cls.calculate_length(child)
            if max_length < length:
                max_length = length
        return max_length


class Tree:
    def __init__(self, root: Node=None):
        self._newick = ''
        self.root = root  # type: Node

    def get_nodes(self) -> List[Node]:
        nodes = []  # type: List[Node]
        queue = [self.root]  # type: List[Node]

        index = 0
        while index < len(queue):
            current = queue[index]
            nodes.append(current)
            if not current.is_leaf():
                queue.extend(current.children)
            index += 1
        return nodes

    def get_internal_nodes(self) -> List[Node]:
        internal_nodes = []  # type: List[Node]
        queue = [self.root]  # type: List[Node]

        index = 0
        while index < len(queue):
            current = queue[index]
            if not current.is_leaf():
                internal_nodes.append(current)
                queue.extend(current.children)
        return internal_nodes

    def get_leaves(self) -> List[Node]:
        leaves = []  # type: List[Node]
        self.bfs(self.root, leaves)
        return leaves

    @classmethod
    def bfs(cls, parent: Node, nodes: List[Node]) -> None:
        if parent.is_leaf():
            nodes.append(parent)
        else:
            for child in parent.children:
                cls.bfs(child, nodes)

    def export_as_newick(self, scaled=True) -> str:
        if self.root is None:
            return ''
        return self.traversal(self.root, scaled)

    @classmethod
    def traversal(cls, node: Node, scaled=True) -> str:
        if node.is_leaf():
            return (node.name + ':' + str(node.distance)) if scaled else node.name

        result = ''
        for i, child in enumerate(node.children):
            result = result + cls.traversal(child)
            if i + 1 < len(node.children):
                result = result + ','

        result = '(' + result + ')'
        label = node.name
        if node.is_root():
            label = label + ';'
        elif scaled:
            label = label + ':' + str(node.distance)

        return result + label

    @classmethod
    def parse_newick(cls, newick: str) -> Optional['Tree']:
        if not newick.endswith(';'):
            return None
        tree = cls()
        tree._newick = newick

        name = ''
        distance = 0.
        index = newick.rfind(')')
        if index + 1 < len(newick) - 1:
            string = newick[index+1: len(newick)-1]
            separator_idx = string.rfind(':')
            if separator_idx < 0:
                name = string
            else:
                name = string[:separator_idx]
                distance = float(string[separator_idx+1:])

        root = Node(name)
        root.distance = distance
        tree.root = root

        parents = [root]  # type: List[Node]
        queue = [newick[1: index]]  # type: List[str]
        indent = '\t'

        q_idx = 0
        while q_idx < len(queue):
            parent = parents[q_idx]
            current = queue[q_idx]
            i = 0
            while i < len(current):
                if current[i] == '(':
                    end = 0
                    j = i + 1
                    while j < len(current):
                        if current[j] == ')':
                            if end != 0:
                                end -= 1
                                j += 1
                                continue

                            subtree = current[i+1:j]
                            if len(current) > j+1:
                                end = current.find(',', j+1)
                                if end < 0:
                                    values = current[j+1:]
                                    i = len(current)
                                else:
                                    values = current[j+1:end]
                                    i = end
                            else:
                                values = ''
                                i = j

                            node_distance = 0.
                            dist_idx = values.rfind(':')
                            if dist_idx < 0:
                                node_name = values
                            else:
                                node_name = values[:dist_idx]
                                try:
                                    node_distance = float(values[dist_idx+1:])
                                except ValueError as e:
                                    from traceback import print_exc
                                    print_exc()

                            node = Node(name=node_name)
                            parent.add_child(node, distance=node_distance)
                            parents.append(node)
                            queue.append(subtree)
                            break
                        elif current[j] == '(':
                            end += 1
                        j += 1
                else:
                    end = current.find(',', i)
                    if end < 0:
                        values = current[i:]
                        node_distance = 0.
                        end = values.rfind(':')
                        if end < 0:
                            node_name = values
                        else:
                            node_name = values[:end]
                            try:
                                node_distance = float(values[end+1:])
                            except ValueError as e:
                                from traceback import print_exc
                                print_exc()
                        node = Node(name=node_name)
                        parent.add_child(node, distance=node_distance)
                        break

                    values = current[i:end]
                    i = end
                    node_distance = 0.
                    end = values.rfind(':')
                    if end < 0:
                        node_name = values
                    else:
                        node_name = values[:end]
                        try:
                            node_distance = float(values[end+1:])
                        except ValueError as e:
                            from traceback import print_exc
                            print_exc()
                    node = Node(name=node_name)
                    parent.add_child(node, distance=node_distance)
                    if end + 1 == len(current):
                        last_node = Node(name='')
                        parent.add_child(last_node, distance=0.)

                i += 1

            indent = indent + '\t'
            q_idx += 1

        return tree

    # @classmethod
    # def import_newicks(cls, newick_multiple: str) -> List['Tree']:
    #     trees = []  # type: List[Tree]
    #
    #     offset = 0
    #     for i, char in newick_multiple:
    #         if char == ';':
    #             string = newick_multiple[offset: i+1].strip()
    #             trees.append(cls.parse_newick(string))
    #             offset = i + 1
    #
    #     return trees

    @classmethod
    def import_newicks(cls, newick_multiple: Union[str, bytes]) -> List['Tree']:
        trees = []  # type: List[Tree]

        offset = 0
        while True:
            try:
                i = newick_multiple.index(';', offset)
            except ValueError:
                break
            string = newick_multiple[offset: i+1].strip()
            trees.append(cls.parse_newick(string))
            offset = i + 1

        return trees

    @classmethod
    def print_dfs(cls, indent: str, parent: Node) -> None:
        print(indent + '"' + parent.name + '"' + ' (' + str(parent.distance) + ')')

        for child in parent.children:
            cls.print_dfs(indent + '  ', child)


class TreeMetaData:
    def __init__(self):
        self.tree = None  # type: Tree
        self.species_check = []  # type: List[bool]
        self.missing_species = []  # type: List[int]


class Species:
    def __init__(self, name: str):
        self.name = name
        self.count = 0
        self.index = 0


class Branch:
    def __init__(self, branch_tag: 'BranchTag'):
        self.branch_tag = branch_tag
        self.index = 0
        self.score = 0.
        self.left_score = 0.
        self.right_score = 0.


class BranchTag:
    LEFT = 0
    RIGHT = 1
    MISSING = 2
    use_hash_array = False

    def __init__(self, size: int, species_on_right: List[int], species_missing: List[int]):
        # produceIdFromIndexedList
        self._id = [0]*size
        for index in species_on_right:
            self._id[index] = 1
        for index in species_missing:
            self._id[index] = 2

        if species_on_right and size != len(species_on_right):
            self._is_root = False
        else:
            self._is_root = True

        if len(species_on_right) != 1 and size - len(species_on_right) != 1:
            self._is_leaf = False
        else:
            self._is_leaf = True

        self._hash_code = 0
        self._tag_length = 0
        self._hash_length = 0
        self._hash_array = []  # type: List[int]
        self._calc_hash_code()

    def get_copy(self, outgroup_index: int=None, side: bool=None):
        copy = BranchTag.__new__(BranchTag)
        copy._id = []
        copy._is_leaf = self._is_leaf
        copy._is_root = self._is_root
        copy._hash_code = self._hash_code

        if side:
            raise NotImplementedError

        if self._is_root or self._id[outgroup_index] == 0:
            copy._id = self._id[:]
        else:
            for x in self._id:
                if x == 0:
                    copy._id.append(1)
                elif x == 1:
                    copy._id.append(0)
                else:
                    copy._id.append(x)

        return copy

    def get_id(self):
        return self._id[:]

    @property
    def is_root(self):
        return self._is_root

    @property
    def is_leaf(self):
        return self._is_leaf

    def get_relationship_with(self, other: 'BranchTag'):
        raise NotImplementedError

    def __copy__(self):
        branch = BranchTag.__new__(BranchTag)
        branch._id = self._id[:]
        branch._is_root = self._is_root
        branch._is_leaf = self._is_leaf
        branch._hash_code = self._hash_code
        return branch

    def __eq__(self, other):
        if not isinstance(other, BranchTag):
            return False

        if self._tag_length != other._tag_length:
            return False
        if self._tag_length == 0 and other._tag_length == 0:
            return True

        if self.use_hash_array:
            for index in range(self._hash_length):
                if self._hash_array[index] != other._hash_array[index]:
                    return False
        else:
            index = 0
            while index < len(self._id) and self._id[index] == 2 and other._id[index] == 2:
                index += 1

            if index < len(self._id):
                first_xor = self._id[index] != other._id[index]
                for i in range(len(self._id)):
                    if self._id[i] != 2 and other._id[i] != 2:
                        if first_xor != (self._id[i] != other._id[i]):
                            return False
                    elif self._id[i] != other._id[i]:
                        return False
        return True

    def equals_ignore_missing(self, obj, missing_tolerance=0, mismatch_tolerance=0):
        if not isinstance(obj, BranchTag):
            return False

        if missing_tolerance < 0:
            missing_tolerance = len(self._id)
        if mismatch_tolerance < 0:
            mismatch_tolerance = len(self._id)

        if len(self._id) != len(obj._id):
            return False
        elif len(self._id) == 0 and len(obj._id) == 0:
            return True

        missing = 0
        mismatch = 0
        index = 0
        while index < len(self._id) and (self._id[index] == 2 or obj._id[index] == 2):
            if self._id[index] != 2 or obj._id[index] != 2:
                missing += 1
                if missing > missing_tolerance:
                    return False
            index += 1

        if index < len(self._id):
            first_xor = self._id[index] != obj._id[index]
            for i in range(index+1, len(self._id)):
                if self._id[i] != 2 and obj._id[i] != 2:
                    if first_xor != (self._id[i] != obj._id[i]):
                        mismatch += 1
                        if mismatch > mismatch_tolerance:
                            return False
                elif self._id[index] != 2 or obj._id[index] != 2:
                    missing += 1
                    if missing > missing_tolerance:
                        return False

        return True

    def equals_side_sensitive(self, obj):
        raise NotImplementedError

    def _calc_hash_code(self):
        from math import ceil
        self._tag_length = len(self._id)
        remainder = 20
        current = 0
        if self.use_hash_array:
            self._hash_length = max(1, int(ceil(self._tag_length/20.)))
            self._hash_array = [int()]*self._hash_length
            remainder -= self._hash_length*20 - self._tag_length

        self._hash_code = 0
        if self._id:
            is_first = True
            first = 2
            for next_ in self._id:
                if next_ == 2:
                    point = 1 if is_first else 2
                else:
                    if is_first:
                        is_first = False
                        first = next_
                    point = 0 if first == next_ else 1

                self._hash_code = 3*self._hash_code + point
                if self.use_hash_array:
                    self._hash_array[current] = 3*self._hash_array[current] + point
                    remainder -= 1
                    if remainder == 0:
                        remainder = 20
                        current += 1

    def __hash__(self):
        return self._hash_code

    def __str__(self):
        return str(self._id)

    def to_consense_string(self):
        raise NotImplementedError

    @classmethod
    def main(cls):
        raise NotImplementedError


class BranchAnalysis:
    def __init__(self, trees):
        self._species_list = []  # type: List[Species]
        self._species_map = {}  # type: Dict[str, Species]
        self._branch_list = []  # type: List[Branch]
        self._branch_map = {}  # type: Dict[BranchTag, Branch]
        self._max_score = float()
        self._trees = self._as_tree_list(trees)  # type: List[Tree]
        self.tree_meta_list = []  # type: List[TreeMetaData]
        self.local_branch_maps = []  # type: List[Dict[BranchTag, List[bool]]]
        self._run_analysis()

    @staticmethod
    def _as_tree_list(trees) -> List[Tree]:
        if hasattr(trees, 'readable') and hasattr(trees, 'read'):  # concatenated newick trees (file)
            return Tree.import_newicks(trees.read())
        elif isinstance(trees, str) or isinstance(trees, bytes):  # concatenated newick trees (string)
            return Tree.import_newicks(trees)
        elif all(isinstance(tree, Tree) for tree in trees):  # [Tree, Tree, ...]
            return trees
        elif all(isinstance(tree, str) for tree in trees):  # [newick, newick, ...]
            tree_list = []
            for s in trees:
                tree = Tree.parse_newick(s)
                if tree:
                    tree_list.append(tree)
                else:
                    raise ValueError
            return tree_list
        else:
            raise ValueError

    def _run_analysis(self) -> None:
        self._max_score = 0.
        self._get_species()
        self._get_branches()

    def _get_species(self) -> None:
        for tree in self._trees:
            leaves = tree.get_leaves()
            for leaf in leaves:
                name = leaf.name
                if name not in self._species_map:
                    species = Species(name)
                    species.count = 1
                    self._species_map[name] = species
                    self._species_list.append(species)
                else:
                    species = self._species_map[name]
                    species.count += 1

        self._species_list.sort(key=lambda sp: (sp.count, sp.name))

        for i, species in enumerate(self._species_list):
            species.index = i

    def _get_branches(self) -> None:
        for tree in self._trees:
            tree_metadata = TreeMetaData()
            local_branch_map = {}  # type: Dict[BranchTag, List[bool]]
            leaves = tree.get_leaves()
            n_leaves = len(leaves)
            total_branches = len(tree.get_nodes())

            existing_species = [False]*len(self._species_list)
            for leaf in leaves:
                existing_species[self._species_map[leaf.name].index] = True
            tree_metadata.species_check = existing_species[:]

            missing_species = []  # type: List[int]
            for j in range(len(existing_species)):
                if not existing_species[j]:
                    missing_species.append(j)
            tree_metadata.missing_species = missing_species[:]

            j = len(self._species_list) - len(leaves)
            points = 1.
            self._max_score += points  # TODO points?

            queue = [tree.root]
            current = 0
            while current < len(queue):
                current_node = queue[current]
                current_leaves = Tree(current_node).get_leaves()
                species_on_right = []  # type: List[int]
                for leaf in current_leaves:
                    species = self._species_map[leaf.name]
                    species_on_right.append(species.index)

                branch_tag = BranchTag(len(self._species_list), species_on_right, missing_species)
                index = 0
                while index < len(branch_tag.get_id()) and branch_tag.get_id()[index] == 2:
                    index += 1

                left_and_right = local_branch_map.get(branch_tag)
                if left_and_right is None:
                    left_and_right = [bool()]*2
                    if index < len(branch_tag.get_id()):
                        if branch_tag.get_id()[index] == 0:
                            left_and_right[0] = True
                        else:
                            left_and_right[1] = True
                    else:
                        left_and_right[0] = True
                        left_and_right[1] = True

                    local_branch_map[branch_tag] = left_and_right
                    branch = self._branch_map.get(branch_tag)
                    if branch is None:
                        branch = Branch(branch_tag)
                        self._branch_map[branch_tag] = branch
                        branch.score = points
                        self._branch_list.append(branch)
                    else:
                        branch.score += points

                    if left_and_right[0]:
                        branch.left_score += points
                    if left_and_right[1]:
                        branch.right_score += points

                else:
                    branch = self._branch_map.get(branch_tag)
                    if index < len(branch_tag.get_id()):
                        if branch_tag.get_id()[index] == 0:
                            if not left_and_right[0]:
                                left_and_right[0] = True
                                branch.left_score += points
                        elif not left_and_right[1]:
                            left_and_right[1] = True
                            branch.right_score += points
                    else:
                        if not left_and_right[0]:
                            left_and_right[0] = True
                            branch.left_score += points
                        if not left_and_right[1]:
                            left_and_right[1] = True
                            branch.right_score += points

                for child in current_node.children:
                    queue.append(child)

                current += 1

            self.local_branch_maps.append(local_branch_map)
            self.tree_meta_list.append(tree_metadata)

        pass
        for i, branch in enumerate(self._branch_list):
            branch.index = i

    @classmethod
    def _propagate(cls, charset: List[str], k: int, position: int, string: str, strings: List[str]) -> None:
        raise NotImplementedError

    @classmethod
    def get_all_possible_kmers(cls, charset, k) -> List[str]:
        raise NotImplementedError

    def print_all_species(self) -> None:
        raise NotImplementedError

    @staticmethod
    def _species_list_to_text_tables(species_list) -> str:
        raise NotImplementedError

    @staticmethod
    def _branch_list_to_text_table(species_list, branch_list, max_score) -> str:
        raise NotImplementedError

    def get_analysis(self) -> str:
        raise NotImplementedError

    def get_branch_information(self, branch_tag, species_on_right, missing_species):
        raise NotImplementedError

    @classmethod
    def _grow_branches(cls, species_list, node):
        raise NotImplementedError

    @staticmethod
    def _build_consensus_tree(max_score, species_list, branch_list):
        raise NotImplementedError

    def get_rooted_branches(self):
        raise NotImplementedError

    def create_consensus_tree(self, outgroup_index=None):
        if outgroup_index is None:
            new_branch_list = self.get_rooted_branches()
        elif outgroup_index >= len(self._species_list):
            return None
        else:
            new_branch_list = []
            for branch in self._branch_list:
                new_branch = Branch(branch.branch_tag.get_copy(outgroup_index))
                new_branch.score = branch.score
                new_branch_list.append(new_branch)

        tree = self._build_consensus_tree(self._max_score, self._species_list, new_branch_list)
        return tree.export_as_newick(scaled=False)

    def mark_tree(self, tree: Union[Tree, str, IO], *,
                  is_rooted: bool, is_core: bool,
                  missing_tolerance: int, mismatch_tolerance: int) -> Optional[str]:
        if isinstance(tree, IO) or hasattr(tree, 'read'):
            tree = Tree.parse_newick(tree.read())
        elif isinstance(tree, str):
            tree = Tree.parse_newick(tree)
            if tree is None:
                return

        leaves = tree.get_leaves()

        existing_species = [bool()]*len(self._species_list)
        for i, leaf in enumerate(leaves):
            if leaf.name in self._species_map:
                j = self._species_map[leaf.name].index
                existing_species[j] = True

        missing_species = []
        for i, exists in enumerate(existing_species):
            if not exists:
                missing_species.append(i)

        queue = [tree.root]
        current = 0
        while current < len(queue):
            current_node = queue[current]
            current_leaves = Tree(current_node).get_leaves()
            is_invalid_branch = False
            species_on_right = []  # type: List[int]

            for leaf in current_leaves:
                if leaf.name not in self._species_map:
                    is_invalid_branch = True
                    if not is_core:
                        current_node.name = '{:.2f}'.format(0)
                    else:
                        current_node.name = '{:d}'.format(0)
                    break
                index = self._species_map[leaf.name].index
                species_on_right.append(index)

            if not is_invalid_branch:
                score = 0.
                for i in range(len(self.tree_meta_list)):
                    tmd = self.tree_meta_list[i]
                    use_hash = True
                    for sp in missing_species:
                        if tmd.species_check[sp]:
                            use_hash = False
                            break

                    if missing_tolerance >= 0 or mismatch_tolerance != 0:
                        use_hash = False

                    if use_hash:
                        branch_tag = BranchTag(len(self._species_list), species_on_right, tmd.missing_species)
                        left_and_right = self.local_branch_maps[i].get(branch_tag)
                        if left_and_right is not None:
                            score += 1.
                    else:
                        branch_tag = BranchTag(len(self._species_list), species_on_right, missing_species)
                        _branch_tags_iterator = self.local_branch_maps[i].keys()
                        for tag in _branch_tags_iterator:
                            if tag.equals_ignore_missing(branch_tag, missing_tolerance, mismatch_tolerance):
                                score += 1.
                                break

                if not current_node.is_leaf():
                    if not is_core:
                        current_node.name = '{:.2f}'.format(score/self._max_score)
                    else:
                        current_node.name = '{:d}'.format(int(score))

            for child in current_node.children:
                queue.append(child)

            current += 1

        return tree.export_as_newick()

    def get_species_list(self):
        return self._species_list[:]

    @staticmethod
    def print_help_message():
        print('Branch Analysis 0.0.1.dev0\n')
        print('Usage:')
        print('    Interactive Mode: \n\t'
              'branchanalysis.py -i')
        print('    Consensus Tree Creation: \n\t'
              'branchanalysis.py -c <input_file> <outgroup_index>')
        print('    Bootstrap Value Marking:	\n\t'
              'branchanalysis.py -b <input_file> <target_newick_file> [-r])')
        print()
        print('<input_file>\t\tName of the file that contains all the newick trees."')
        print('<outgroup_index>\tIndex of the item to be used as the outgroup.')
        print('\t\t\t\tIf this value is 0, there will be no outgroup \n'
              '\t\t\t\tand the input trees will be treated as rooted trees.')
        print('<target_newick_file>\tThe newick file to be marked.')
        print('-r\t\t\tUse if the target tree is rooted.')

    @classmethod
    def run_test(cls):
        raise NotImplementedError

    @classmethod
    def run_interactive_mode(cls):
        raise NotImplementedError

    @classmethod
    def main(cls, *args):
        if len(args) == 0:
            cls.print_help_message()
        elif args[0] == '-i':
            print('!!! Not Implemented !!!', file=sys.stderr)
            cls.run_interactive_mode()
        else:
            if args[0] == '-c':
                input_file = open(args[1])
                branch_analysis = cls(input_file)
                outgroup_index = int(args[2])
                if outgroup_index == 0:
                    marked_tree = branch_analysis.create_consensus_tree()
                else:
                    marked_tree = branch_analysis.create_consensus_tree(outgroup_index)
                print('Consensus tree created', file=sys.stderr)
                print(marked_tree)
            elif args[0] == '-b':
                input_file = open(args[1])
                branch_analysis = cls(input_file)
                target_newick_file = open(args[2])
                if len(args) == 4 and args[3] == '-r':
                    marked_tree = branch_analysis.mark_tree(target_newick_file,
                                                            is_rooted=True, is_core=False,
                                                            missing_tolerance=-1, mismatch_tolerance=0)
                elif len(args) == 4 and args[3] == '-bcg':
                    marked_tree = branch_analysis.mark_tree(target_newick_file,
                                                            is_rooted=False, is_core=True,
                                                            missing_tolerance=-1, mismatch_tolerance=0)
                else:
                    marked_tree = branch_analysis.mark_tree(target_newick_file,
                                                            is_rooted=False, is_core=False,
                                                            missing_tolerance=-1, mismatch_tolerance=0)
                print(marked_tree)
            elif args[0] == '-t':
                cls.run_test()
            else:
                cls.print_help_message()


if __name__ == '__main__':
    BranchAnalysis.main(sys.argv[1:])
