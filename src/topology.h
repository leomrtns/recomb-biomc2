/*! 
 * Copyright (C) 2006	Leonardo de Oliveira Martins
 * 
 * leo at lbm ab a u-tokyo ac jp 
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the 
 * Free Software Foundation, Inc., 51 Franklin Street, 
 * Fifth Floor, Boston, MA  02110-1301, USA.
 */

/*! \defgroup topology Topology functions/structures */

/*! \file topology.h 
 *  \ingroup topology
 *  \brief Header file for topology.c, topology_distance.c and topology_split_partition.c .
 *  
 */

#ifndef _biomc2_topology_h_
#define _biomc2_topology_h_

#include "nexus.h"

typedef struct topol_node_struct* topol_node;
typedef struct topology_struct* topology;
typedef struct split_space_struct* split_space;
typedef struct split_partition_struct split_partition;

/*! \brief Compare two topologies based on their splits, here represented as
 * bit-strings (e.g. the pair 00110001, 11001110). */
struct split_space_struct 
{
	/*! \brief Number of internal nodes. */
	int inodes;
	/*! \brief Number of leaves. */
	int n_leaves;
	/*! \brief Vector with pointers to split_partitions of topology 1. */
	split_partition **u1; 
	/*! \brief Vector with pointers to split_partitions of topology 2. */
	split_partition **u2;
	/*! \brief Split_partition vector size of topology 1. */
	int n_u1; 
	/*! \brief Split_partition vector size of topology 2. */
	int n_u2;
	/*! \brief Index of original leaves mapped to compressed leaves */
	int *leaf_idx;
	/*! \brief Leaves involved in recombination */
	int *leaf_recomb;
	/*! \brief Number of active leaves. */
	int n_active;
	/*! \brief Edges in agreement. */
	split_partition **agree;
	/*! \brief Number of edges in agreement. */
	int n_agree;
	/*! \brief List of possible prune subtrees in one direction. */
	split_partition *dis1;
	/*! \brief List of possible prune subtrees in oposite direction. */
	split_partition *dis2;
	/*! \brief List of pointers to possible prune subtrees (smallest between dis1 and dis2). */
  split_partition **d_ptr;
	/*! \brief Number of possible prune subtrees. */
	int n_disagree;
	/*! \brief Smallest XOR in agreement prune subtree. */
	split_partition *best_prune;
	/*! \brief Active leaves in mask notation (e.g. 001111). */
  long long int active_l;
	/*! \brief SPR distance, and RF distance */
	int dSPR, rfdistance;
	/*! \brief Pointer to opaque function that stores pruned subtrees (or does nothing) */ 
	void (*extend_subtree_info) (int *index, int *idx_n, split_space split);
	/*! \brief Pointer to opaque function that maps original leaves to compressed leaves (or does nothing) */
	void (*store_compression_info) (int *to, int *from, split_space split);
};

/*! \brief Bit-string representation of splits. */
struct split_partition_struct 
{
	/*! \brief Boolean vector with bipartition (poor man's bit string) (unused). */
	int b[MaxLeaves];
	/*! \brief Integer representation of bipartition (bit string) where "1" is the
	 * side with fewer leaves (smaller subtree). */
  long long int l;
  /*! \brief Number of leaves present in smaller bipartition. */
	int n;
	/*! \brief Back-reference to node. */
	topol_node node;
};

/*! \brief Information of a node (binary tree). */
struct topol_node_struct
{
	/*! \brief Parent node. */
	topol_node up;
	/*! \brief Child node. */
	topol_node right;
	/*! \brief Child node. */
	topol_node left;
	/*! \brief Sister node. */
	topol_node sister;
	/*! \brief Node ID (values smaller than nleaves indicate leaves). */
	int id;
	/*! \brief Mapping between nodes of different topologies, used by swap_likelihood(). */
	int map_id;
	/*! \brief If internal node, TRUE; if leaf, FALSE. */ 
	bool internal;
	/*! \brief Quartet partitioning (bipartition) induced by edge (leaf set rooted at node). */
	split_partition bipartition;
	/*! \brief Are partial likelihoods (upper part) already calculated? */
	bool u_done; 
	/*! \brief Are partial likelihoods (lower part) already calculated? */
	bool d_done;
};

/*! \brief Binary unrooted topology (rooted at leaf with ID zero) */
struct topology_struct
{
	/*! \brief Topology ID (in linked_list_struct), used only for debug purposes. */
	int id;
	/*! \brief Vector with pointers to nodes. */
	topol_node *nodelist;
	/*! \brief Pointer to root node (a leaf, since the tree in unrooted). */
	topol_node root;
	/*! \brief Number of leaves. */
	int nleaves;
	/*! \brief Number of nodes, including leaves. Since the tree is binary and
	 * unrooted, the number of nodes equals \f$ 2L-2\f$ where \f$L\f$ is the 
	 * number of leaves. */
	int nnodes;
	/*! \brief How to revert SPR move (prune node). */
	topol_node undo_prune; 
	/*! \brief How to revert SPR move (regraft node). */
	topol_node undo_regraft;

	/*! \brief Pointer to previous topology in linked list. */
 	topology prev;
	/*! \brief Pointer to next topology in linked list. */
	topology next;
	/*! \brief First segment (over the alignment) for this topology. */
	int first_segment; 
	/*! \brief Last segment (over the alignment) for this topology. */
	int last_segment; 
	
	double likelihood_accepted, /*!< \brief Accepted likelihood of region (e.g. after an iteration). */
				 likelihood_current,  /*!< \brief Current likelihood of region (e.g. inside Al-Awadhi update, 
																it is only temporarily accepted). */
				 likelihood_proposal; /*!< \brief Proposal likelihood of region (calculated whenever a 
																parameter/topology structure changes). */

	/*! \brief Posterior frequency of topology (used by OLD VERSION OF summarise.c to find MAP topology). */
	int tfreq;
};

/* Global function prototypes */

/* topology.c */

/*! \brief Allocate space for new topology_struct */
topology new_topology (int nleaves);

/*! \brief Free space allocated by topology_struct */
void del_topology (topology topol);

/*! \brief Copy information from nexus_tree struct to topology struct 
 *
 * Since topology nodes are related to likelihood vectors their IDs follow
 * strict rules: 
 * - their IDs should not change, only the relations between them represented by
 *   pointers up, left, right, sister;
 * - their IDs are the position in the topology_struct::nodelist vector (with actual nodes);
 * - IDs smaller than number of leaves represent the leaves (indexed by same
 * IDs);
 * \param[in] nxs_tree nexus_tree_struct with node IDs respecting topology rules
 * \param[out] tree (previously allocated) copied topology_struct
 */
void copy_topology_from_nexus_tree (topology tree, nexus_tree nxs_tree);

/*! \brief Copy information from topology_struct, preserving possible common
 * subtrees. 
 * 
 * Calls copy_topology_from_topology() after mapping topol_node_struct between
 * original and copied topologies; The mapping is accomplished by
 * - updating topol_node_struct::map_id if the split or the ID are the same
 * - setting topol_node_struct::u_done or topol_node_struct::d_done to TRUE if
 * corresponding subtrees are the same. 
 * This can speed up likelihood calculations when both topologies have subtrees
 * in common - in which case swap_likelihood() must be called to invert
 * corresponding nodes.
 * \param[in]  split split_space_struct necessary to compare topologies splits
 * \param[in]  from_tree original topology_struct 
 * \param[out] to_tree (previously allocated) copied topology_struct
 */
void copy_topology_from_topology_mapping (topology to_tree, topology from_tree, split_space split);

/*! \brief Copy information from topology_struct. 
 *
 * Since IDs do not change, this function only needs to update topol_node_struct::up,
 * topol_node_struct::right, and topol_node_struct::left pointers and topol_node_struct::map_id from
 * internal nodes; update of topol_node::sister is handled by function
 * update_topology_sisters(). 
 * \param[in]  from_tree original topology_struct 
 * \param[out] to_tree (previously allocated) copied topology_struct
 */
void copy_topology_from_topology (topology to_tree, topology from_tree);

/*! \brief Update pointers to topol_node_struct::sister */
void update_topology_sisters (topology tree);

/*! \brief Boolean if node2 is on the path of node1 to the root. */
bool node1_is_child_of_node2 (topol_node node1, topol_node node2);

/*! \brief Reroot topology_struct on new leaf. This function is not used and may be broken */
topol_node reroot_topology (topology tree, topol_node newroot);

/*! \brief Print subtree in newick format to string using leaf IDs.
 *
 * Stores in string the tree in newick format, using leaf ID numbers (needs a
 * TRANLATION nexus block).
 * Memory allocation is handled by this function, but needs to be freed by the calling 
 * function. 
 * \param[in] tree tree to be printed
 * \return a pointer to newly allocated string 
 */
char * topology_to_string (const topology tree);

/*! \brief Print subtree in newick format to string using leaf names.
 *
 * Stores in string the tree in newick format, preserving sequence names.
 * Memory allocation is handled by this function, but needs to be freed by the calling 
 * function. 
 * \param[in] tree tree to be printed
 * \param[in] taxlabel vector of leaf names
 * \return a pointer to newly allocated string 
 */
char * topology_to_string_with_name (const topology T, char **taxlabel);

/*! \brief Print subtree in newick format to string using leaf names with
 * random branch lengths.
 *
 * Stores in string the tree in newick format, generating branch lengths from a
 * uniform between MIN_LENGTH and (MAX_LENGTH + MIN_LENGTH). 
 * Memory allocation is handled by this function, but needs to be freed by the calling 
 * function. 
 * \param[in] tree tree to be printed
 * \param[in] taxlabel vector of leaf names
 * \param[in] min_length minimum branch length
 * \param[in] max_length (maximum - minimum) branch length
 * \return a pointer to newly allocated string 
 */
char * topology_to_string_with_name_generate_length (const topology T, char **taxlabel, 
																										 double min_length, double max_length);


/* topology_distance.c */

/*! \brief Post-roder traversal where bipartition information in 
 * topol_node_struct::split_partition is updated. */
void create_topology_bipartition_l (topology tree);

/*! \brief Store bipartition information in post-order traversal */ 
void postorder_mapping_bipartition (topol_node this, split_partition **u, int *n);

/*! \brief Rearrange topol_node_struct::left and topol_node_struct::right where
 * smaller leaf is always on the left. */
void order_topology_by_bipartition (topol_node this);

/*! \brief Compare topologies and updates equal subtrees. Used by copy_topology_from_topology_mapping(). */
void update_topology_from_topology_split (topology to_tree, topology from_tree, split_space split);

/*! \brief fast comparison between two trees to check if they are the same or not */
bool dSPR_is_zero_l (topology t1, topology t2, split_space split);

/*! \brief Calculate approximate \f$d_{SPR}\f$ between topologies t1 and t2. */
int dSPR_topology_l (topology t1, topology t2, split_space split, int *max);

/*! \brief Copy topol_node_struct::id to vector in pre-order traversal. */
void euler_semitour_mapping (topol_node this, int *u, int *n);

/* topology_split_partition.c */ 

/*! \brief Allocate space for new split_space_struct. */
split_space new_split_space (int n_leaves);

/*!\brief Free memory allocated by split_space_struct. */
void del_split_space (split_space split);

/*! \brief Reset values of split_space_struct before new split bipartition calculation. */
void initialize_split_space_l (split_space split);

/*! \brief Extends \f$ d_{SPR} \f$ to give also the pruned leaves. Should be
 * called before each dSPR_topology_l() since it also initializes the vectors. */
void setup_split_store_recomb_info (split_space split);

/*! \brief \f$ d_{SPR} \f$ extension is not necessary. Can be called just once
 * (and is called by new_split_space() ). */ 
void setup_split_dont_store_recomb_info (split_space split);

/*! \brief Set bipartition of a leaf. */ 
void split_initialize_l (split_partition *bpt, int label_id);

/*! \brief Set bipartition of internal node. */
void split_OR_l (split_partition *res, const split_partition *b1, const split_partition *b2);

/*! \brief Set up split_partition_struct as being the smallest subtree it represents (at first iteration). */
void split_minimize_initial_subtree_size_l (split_partition **u, int *n, split_space split);

/*! \brief Set up split_partition_struct as being the smallest subtree it represents (general use). */
void split_minimize_subtree_size_l (split_partition **u, int *n, split_space split);

/*! \brief Search for common bipartitions, storing them in split_space_struct::agree. */
void split_create_agreement_list_l (split_space split);

/*! \brief Search for common bipartitions, storing them in split_space_struct::agree and 
 * updating topol_node_struct::map_id. */
void split_create_agreement_list_mapping (split_space split);

/*! \brief Remove edges (bipartitions) that are present in split_space_struct:agree. */
void split_remove_agreement_edges_l (split_partition **u, int *n, split_space split);

/*! \brief Remove redundant leaves in bipartitions when there is a subtree of
 * size 2, updating topol_node_struct::u_done or topol_node_struct::d_done. */
void split_compress_subtree_agreement_topol (split_space split);

/*! \brief Remove redundant leaves in bipartitions when there is a subtree of size 2. */
void split_compress_subtree_agreement_l (split_space split);

/*! \brief Calculate disagreement between all bipartition pairs, unless a
 * disagreement of size one (one leaf) is found. */
void split_create_disagreement_list_l (split_space split);

/*! \brief reduce split_space_struct::d_ptr size (created by split_create_disagreement_list_l())
 * by removing identical bipartitions. */
void split_remove_duplicate_disagreement_l (split_space split);

/*! \brief Compare split_space_struct::d_ptr (disagreement) and split_space_struct::agree 
 * (common subtrees) and search for best match. */
void split_find_smallest_disagreement_subtree_l (split_space split);

/*! \brief Reduce disagreement edges subtrees by removing partial subtrees in
 * agreement (unused). */
void split_minimize_disagreement_l (split_space split);

/*! \brief Remove leaves belonging to pruned subtree (smallest subtree in
 * disagreemtent. */
void split_remove_pruned_subtree_l (split_partition *prune, split_space split);

/*! \brief Print to stdout bit-string representation of bipartitition */
void split_print_binary_bipartition_l (const split_partition *u, split_space split);

/*! \brief Count number of leaves in split_partition_struct::l. */
void split_bipartition_size_l (int *count, const long long int  *l, split_space split);

#endif
