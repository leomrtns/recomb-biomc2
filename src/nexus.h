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

/*! \defgroup nexus Nexus file handling functions/structures. */

/*! \file nexus.h 
 *  \ingroup nexus
 * \brief Header file for nexus.c, nexus_tree.c and nexus_alignment.c  
 */

#ifndef _biomc2_nexus_h_
#define _biomc2_nexus_h_

#include "hashtable.h"

/*! \brief maximum name length for taxa (alignment and tree files). */
#define MAX_NAME_LENGTH 128

typedef struct nexus_alignment_struct* nexus_alignment;
typedef struct nexus_node_struct* nexus_node;
typedef struct nexus_tree_struct* nexus_tree;
typedef struct nexus_treespace_struct* nexus_treespace;

/*! \brief Data from alignment file. */
struct nexus_alignment_struct
{
	/*! \brief Number of species, according to sequence file. Tree should not disagree. */
	int ntax;
	/*! \brief Sequence size, assuming it's the same for all taxa. */
	int nchar;
	/*! \brief Vector with aligned sequence for each taxon. */
	char **character;
	/*! \brief Taxon names from file. */
	char **taxlabel;
	/*! \brief Lookup table with taxon names. */
	hashtable taxlabel_hash;
	/*! \brief Number of gene segments (ASSUMPTIONS BLOCK). */
	int n_charset;
	/*! \brief Start of each gene segment (from 1...NCHAR) (ASSUMPTIONS BLOCK). */
	int *charset_start;
	/*! \brief End of each gene segment (ASSUMPTIONS BLOCK). */
	int *charset_end;
};

/*! \brief Node information for each tree in tree file. */
struct nexus_node_struct
{
	/*! \brief Parent node. */
	nexus_node up;
	/*! \brief Child node. */
	nexus_node right;
	/*! \brief Child node. */
	nexus_node left;
	/*! \brief Initial pre-order numbering of node. */
	int id;       
	/*! \brief Branch length from node to node->up. */
	double branch_length;
	/*! \brief Leaf sequence name. Must match sequence name from nexus_alignment_struct. */
	char *taxlabel;
};

/*! \brief Data from each tree in tree file. */
struct nexus_tree_struct
{
	/*! \brief Vector with pointers to every internal node. */
	nexus_node *nodelist;
	/*! \brief Vector with pointers to tree leaves. */
	nexus_node *leaflist;
	/*! \brief Pointer to root node (usually a leaf, since the tree in unrooted). */
	nexus_node root;
	/*! \brief Boolean saying if tree has branch lengths or not. */
	bool has_branches;
	/*! \brief Number of leaves (should correspond to number of sequences in nexus_alignment_struct). */
	int nleaves;
	/*! \brief Number of nodes, including leaves. Since the tree is binary and unrooted, the number of nodes equals 
	 * \f$ 2L-2\f$ where \f$L\f$ is the number of leaves. */
	int nnodes;
};

/*! \brief Collection of trees from tree file. */
struct nexus_treespace_struct
{
	/*! \brief Number of trees in nexus file. */
	int ntrees;
	/*! \brief Vector of nexus_tree. */
	nexus_tree *T;
	/*! \brief Number of leaves (should correspond to number of sequences in nexus_alignment_struct). */
	int nleaves;
	/*! \brief Number of nodes, including leaves. Since the tree is binary and unrooted, the number of nodes equals 
	 * \f$ 2L-2\f$ where \f$L\f$ is the number of leaves. */
	int nnodes;
	/*! \brief Taxon names. Should correspond to sequence names in nexus_alignment_struct. */
	char **taxlabel;
	/*! \brief Lookup table with taxon names. */
	hashtable taxlabel_hash;
};

/* Global function prototypes */

/* nexus.c */
/*! \brief Removes (possible nested/multiline) nexus comments of the form [] (brackets). */
char* remove_nexus_comments (char **string, size_t *stringsize, FILE *stream);
/*! \brief Changes uppercase characters by lowercase versions. */
char* lowercase_string (char *string);
/*! \brief Changes lowercase characters by uppercase versions. */
char* uppercase_string (char *string);
/*! \brief Removes spaces, tabs from string. */
char* remove_space_from_string (char *string);

/* nexus_alignment.c */
/*! \brief Reads DNA alignment from file and store info in nexus_alignment_struct. */
nexus_alignment read_nexus_alignment_file (char *seqfilename);
/*! \brief Prints alignment to screen in FASTA format (debug purposes). */
void print_nexus_alignment (nexus_alignment align);
/*! \brief Frees memory from nexus_alignment_struct. */
void del_nexus_alignment (nexus_alignment align);

/* nexus_tree.c */

/*! \brief Frees memory from nexus_treespace_struct. */
void del_nexus_treespace (nexus_treespace tsp);
/*! \brief Reads tree file and store info in nexus_treespace_struct. */
nexus_treespace read_nexus_treespace_file (char *seqfilename);
/*! \brief Calculates average branch length of nexus_tree_struct. */
double average_rate_nexus_tree (nexus_tree T);
/*! \brief Mapping between align->hashtable and treespace->hashtable. 
 * 
 * Given a vector with names (taxlabel) and a hashtable (hash), we want to number 
 * taxlabel according to hash values. 
 * After checking if all names from taxlabel have a corresponding hash key,
 * this function will create a vector with the position, in hash, of elements of
 * taxlabel. This is what we call a mapping. Ex:
 * If we have the hashtable and taxlabels are
 * \verbatim
 *     hash["A"] = 0      taxlabel[0] = "C" 
 *     hash["B"] = 1      taxlabel[1] = "B" 
 *     hash["C"] = 2      taxlabel[2] = "D" 
 *     hash["D"] = 3      taxlabel[3] = "E" 
 *     hash["E"] = 4      taxlabel[4] = "A" 
 * \endverbatim
 * then we have that the mapping is given by
 * \verbatim
 *     order[0] = 2
 *     order[1] = 1
 *     order[2] = 3
 *     order[3] = 4
 *     order[4] = 0
 * \endverbatim
 * This is necessary since despite nexus_tree_struct has freedom about the order
 * of nodes (including leaves), in topology_struct the order is defined.
 * \param[in]  hash Hashtable from nexus_alignment_struct. 
 * \param[in]  taxlabel Name vector from nexus_tree_struct. 
 * \param[in]  size Number of elements in taxlabel and hash.
 * \param[out] order Mapping between hash and taxlabel.
 */
void map_hashtable_order (int *order, hashtable hash, char **taxlabel, int size);

/*! \brief Relabels tree nodes using mapping ( id = order[id] ).
 *
 * Using the order defined by map_hashtable_order(), all trees belonging to
 * a nexus_treespace_struct will be relabeled and rerooted by the mapping, where
 * the new root will be the leaf having order[]=0. Notice than even in the same
 * nexus_treespace_struct distinct trees may have distinct roots and leaf
 * label orders, despite the taxlabel[] vector with leaf names is shared.
 * \param[in] order Mapping to be applied to nexus_tree_struct elements.
 * \param[in,out] t nexus_treespace_struct with nexus_tree_struct elements to be
 * reordered.
 * */
void order_nexus_treespace_id (nexus_treespace t, int *order);
/*! \brief Reroots nexus_tree_struct on new leaf (so that all trees have same root leaf). */
nexus_node reroot_nexus_tree (nexus_tree tree, nexus_node newroot);
/*! \brief Prints subtree in newick format to string.
 *
 * Stores in string the tree in newick format. Memory allocation is handled by this
 * function, but needs to be freed by the calling function.
 * This is not the best alternative, since there may be a limit on the maximum
 * string size handled by printf() family functions. In some systems this limit
 * is 519 bytes (according to the glibc documentation).
 * \param[in] tree Tree to be printed.
 * \return pointer to string with tree in newick format.
 */
char* nexus_tree_to_string (const nexus_tree tree);
/*! \brief Prints subtree in dot format to file.
 *
 * Prints to file the tree in dot format (undirected graph). 
 * The dot format can be used with the <a href="http://www.graphviz.org/">graphviz</a> 
 * suite of programs, and is not restricted to trees but can also handle
 * arbitrary graph structures. Notice that we do not make use of the graphviz
 * library, we simply create the text file graphviz programs take as input.
 * Unfortunately, it is not helpful to print the nexus_tree_struct since the
 * program works basically with the topology_struct. On the other hand it is
 * easy to change this function to make it work with topology_struct.
 * \param[in,out] fout pointer to file handler where tree is to be printed;
 * \param[in] label graph name or any other label;
 * \param[in] tree tree to be printed;
 */
void graphviz_file_nexus_tree (FILE * fout, char *label, const nexus_tree tree);

#endif
