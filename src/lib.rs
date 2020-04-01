use std::cmp;

use nalgebra::DMatrix;

pub struct Tree<'a> {
    post_order: Vec<&'a TreeNode>,
    left_most_leaf_descendant: Vec<usize>,
    key_roots: Vec<usize>,
}

pub struct TreeNode {
    label: String,
    children: Vec<Box<TreeNode>>,
}

impl TreeNode {
    pub fn new(label: &str) -> TreeNode {
        TreeNode {
            label: String::from(label),
            children: vec![],
        }
    }

    pub fn with_children(mut self, children : Vec<Box<TreeNode>>) -> TreeNode {
        self.children = children;
        self
    }

    fn post_order(&self) -> Vec<&TreeNode> {
        let mut result = vec![];
        for child in self.children.iter() {
            result.extend(child.post_order().iter());
        }
        result.push(self);
        result
    }
}

impl<'a> Tree<'a> {
    pub fn new(root: &'a TreeNode) -> Tree {
        let post_order = root.post_order();
        let left_most_leaf_descendant = Tree::left_most_leaf_descendant(&post_order);
        let key_roots = Tree::keyroots(&post_order);
        Tree {
            post_order,
            left_most_leaf_descendant,
            key_roots,
        }
    }

    fn left_most_leaf_descendant(post_order: &[&TreeNode]) -> Vec<usize> {
        let mut result = Vec::with_capacity(post_order.len());
        for (idx, node) in post_order.iter().enumerate() {
            let left_most_child = idx - (node.post_order().len() - 1);
            result.push(left_most_child);
        }
        result
    }

    fn keyroots(post_order: &[&TreeNode]) -> Vec<usize> {
        let mut key_roots = Vec::new();
        let mut to_look_at = Vec::new();
        // insert root node
        key_roots.push(post_order.len() - 1);
        // insert appropriate child nodes
        to_look_at.push(post_order.len() - 1);
        while !to_look_at.is_empty() {
            let n = to_look_at.pop().unwrap();
            for (idx, child) in post_order[n].children.iter().enumerate() {
                if idx > 0 {
                    key_roots.push(Tree::id(post_order, child));
                }
                to_look_at.push(Tree::id(post_order, child));
            }
        }
        // sort keyroots in ascending order
        key_roots.sort();
        key_roots
    }

    fn id(post_order: &[&TreeNode], node: &TreeNode) -> usize {
        post_order
            .iter()
            .position(|&n| n as *const TreeNode == node as *const TreeNode)
            .unwrap()
    }

    fn label_cmp(tn1: &TreeNode, tn2: &TreeNode, relabeling_cost: u64) -> u64 {
        if tn1.label == tn2.label {
            0u64
        } else {
            relabeling_cost
        }
    }

    fn forest_distance(
        key_root_1: usize,
        key_root_2: usize,
        t1: &Tree,
        t2: &Tree,
        td: &mut DMatrix<u64>,
        insertion_cost: u64,
        deletion_cost: u64,
        relabeling_cost: u64,
    ) {
        let l1_i = t1.left_most_leaf_descendant[key_root_1];
        let l2_j = t2.left_most_leaf_descendant[key_root_2];
        let mut fd: DMatrix<u64> = DMatrix::zeros(key_root_1 - l1_i + 2, key_root_2 - l2_j + 2);
        for i in 1..(key_root_1 - l1_i + 2) {
            fd[(i, 0)] = fd[(i - 1, 0)] + deletion_cost;
        }
        for i in 1..(key_root_2 - l2_j + 2) {
            fd[(0, i)] = fd[(0, i - 1)] + insertion_cost;
        }
        for i in 1..(key_root_1 - l1_i + 2) {
            for j in 1..(key_root_2 - l2_j + 2) {
                // check if t1 and t2 are both trees
                if t1.left_most_leaf_descendant[i + l1_i - 1] == l1_i
                    && t2.left_most_leaf_descendant[j + l2_j - 1] == l2_j
                {
                    fd[(i, j)] = cmp::min(
                        cmp::min(
                            fd[(i - 1, j)] + deletion_cost,
                            fd[(i, j - 1)] + insertion_cost,
                        ),
                        fd[(i - 1, j - 1)]
                            + Tree::label_cmp(
                                t1.post_order[i + l1_i - 1],
                                t2.post_order[j + l2_j - 1],
                                relabeling_cost,
                            ),
                    );
                    td[(i + l1_i - 1, j + l2_j - 1)] = fd[(i, j)];
                }
                // in this case at least t1 or t2 is a forest
                else {
                    fd[(i, j)] = cmp::min(
                        cmp::min(
                            fd[(i - 1, j)] + deletion_cost,
                            fd[(i, j - 1)] + insertion_cost,
                        ),
                        fd[(
                            t1.left_most_leaf_descendant[i + l1_i - 1] - l1_i,
                            t2.left_most_leaf_descendant[j + l2_j - 1] - l2_j,
                        )] + td[(i + l1_i - 1, j + l2_j - 1)],
                    );
                }
            }
        }
    }

    pub fn weighted_tree_edit_distance(
        &self,
        other: &Tree,
        insertion_cost: u64,
        deletion_cost: u64,
        relabeling_cost: u64,
    ) -> u64 {
        let mut td: DMatrix<u64> = DMatrix::zeros(self.post_order.len(), other.post_order.len());
        for &x in self.key_roots.iter() {
            for &y in other.key_roots.iter() {
                Tree::forest_distance(
                    x,
                    y,
                    self,
                    other,
                    &mut td,
                    insertion_cost,
                    deletion_cost,
                    relabeling_cost,
                );
            }
        }
        td[(self.post_order.len() - 1, other.post_order.len() - 1)]
    }

    pub fn tree_edit_distance(&self, other: &Tree) -> u64 {
        self.weighted_tree_edit_distance(other, 1, 1, 1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // macro for more convenient creation of trees
    macro_rules! tree {
        ($r:expr) => {
            TreeNode::new($r)
        };
        ($r:expr,[ $( $c:expr ),* ] )=> {{
            let mut root = TreeNode::new($r);
            $(
                root.children.push(Box::new($c));
            )*
            root
        }};
    }

    #[test]
    fn test_post_order() {
        // example tree given by root_node with node number in post order depicted next to label
        //                             A  5
        //                             +
        //                             |
        //                       +-----+-----+
        //                       |     |     |
        //                       v     v     v
        //                       B 0   C 3   D 4
        //                             +
        //                             |
        //                          +--+--+
        //                          |     |
        //                          v     v
        //                          E 1   F 2

        let root_node = tree!(
            "A",
            [tree!("B"), tree!("C", [tree!("E"), tree!("F")]), tree!("D")]
        );
        assert_eq!("B", root_node.post_order()[0].label);
        assert_eq!("E", root_node.post_order()[1].label);
        assert_eq!("F", root_node.post_order()[2].label);
        assert_eq!("C", root_node.post_order()[3].label);
        assert_eq!("D", root_node.post_order()[4].label);
        assert_eq!("A", root_node.post_order()[5].label);
    }

    #[test]
    fn test_leftmost_leaf_descendant() {
        // example tree given by root_node with node number in post order depicted next to label
        //                             A  5
        //                             +
        //                             |
        //                       +-----+-----+
        //                       |     |     |
        //                       v     v     v
        //                       B 0   C 3   D 4
        //                             +
        //                             |
        //                          +--+--+
        //                          |     |
        //                          v     v
        //                          E 1   F 2

        let root_node = tree!(
            "A",
            [tree!("B"), tree!("C", [tree!("E"), tree!("F")]), tree!("D")]
        );
        let tree = Tree::new(&root_node);

        assert_eq!(0, tree.left_most_leaf_descendant[5]);
        assert_eq!(4, tree.left_most_leaf_descendant[4]);
        assert_eq!(1, tree.left_most_leaf_descendant[3]);
        assert_eq!(2, tree.left_most_leaf_descendant[2]);
        assert_eq!(1, tree.left_most_leaf_descendant[1]);
        assert_eq!(0, tree.left_most_leaf_descendant[0]);
    }

    #[test]
    fn test_key_roots() {
        // example tree given by root_node with node number in post order depicted next to label
        //                             A  5
        //                             +
        //                             |
        //                       +-----+-----+
        //                       |     |     |
        //                       v     v     v
        //                       B 0   C 3   D 4
        //                             +
        //                             |
        //                          +--+--+
        //                          |     |
        //                          v     v
        //                          E 1   F 2

        let root_node = tree!(
            "A",
            [tree!("B"), tree!("C", [tree!("E"), tree!("F")]), tree!("D")]
        );
        let tree = Tree::new(&root_node);

        assert_eq!(2, tree.key_roots[0]);
        assert_eq!(3, tree.key_roots[1]);
        assert_eq!(4, tree.key_roots[2]);
        assert_eq!(5, tree.key_roots[3]);
    }

    #[test]
    fn test_self_distance_is_zero() {
        let tree_1_root_node = tree!(
            "A",
            [
                tree!("B"),
                tree!("C", [tree!("C1"), tree!("C2")]),
                tree!("D")
            ]
        );
        let tree_2_root_node = tree!("X");

        let tree_1 = Tree::new(&tree_1_root_node);
        let tree_2 = Tree::new(&tree_2_root_node);

        // distance between a tree and itself should always be zero
        assert_eq!(0, tree_1.tree_edit_distance(&tree_1));
        assert_eq!(0, tree_2.tree_edit_distance(&tree_2));

        // distance to any tree that is different must not be zero
        assert_ne!(0, tree_1.tree_edit_distance(&tree_2));
        assert_ne!(0, tree_2.tree_edit_distance(&tree_1));
    }

    #[test]
    fn test_distance_with_single_node_trees() {
        let tree_1_root_node = tree!("A");
        let tree_2_root_node = tree!("B");

        let tree_1 = Tree::new(&tree_1_root_node);
        let tree_2 = Tree::new(&tree_2_root_node);

        assert_eq!(1, tree_1.tree_edit_distance(&tree_2));
        assert_eq!(1, tree_2.tree_edit_distance(&tree_1));
    }

    #[test]
    fn test_distance_with_trees() {
        let tree_1_root_node = tree!("A", [tree!("B"), tree!("C"), tree!("D", [tree!("E")])]);
        let tree_2_root_node = tree!("X", [tree!("C"), tree!("Y", [tree!("Z")])]);

        let tree_1 = Tree::new(&tree_1_root_node);
        let tree_2 = Tree::new(&tree_2_root_node);

        assert_eq!(4, tree_1.tree_edit_distance(&tree_2));
        assert_eq!(4, tree_2.tree_edit_distance(&tree_1));
    }

    #[test]
    fn test_weighted_distance() {
        let tree_1_root_node = tree!("A");
        let tree_2_root_node = tree!("B");

        let tree_1 = Tree::new(&tree_1_root_node);
        let tree_2 = Tree::new(&tree_2_root_node);

        assert_eq!(2, tree_1.weighted_tree_edit_distance(&tree_2, 1, 1, 3));
        assert_eq!(2, tree_2.weighted_tree_edit_distance(&tree_1, 1, 1, 3));
    }
}
