/* Copyright (C) 2006  Leonardo de Oliveira Martins
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

/*! \file topology.c
 *  \ingroup topology 
 *  \brief Topology information used by phylogeny_struct (without likelihood etc.) 
 */

#include "topology.h"

/* local funtion prototypes */

/*! \brief Auxiliary function to topology_to_string(). */
void topology_subtree_to_string (char *str, const topol_node node);
/*! \brief Auxiliary function to topology_to_string_with_name(). */
void topology_subtree_to_string_with_name (char *str, const topol_node node, char **taxlabel);
/*! \brief Auxiliary function to topology_to_string_with_name_generate_length(). */
void topology_subtree_to_string_with_name_generate_length (char *str, const topol_node node, char **taxlabel,
																													 double *min_length, double *max_length);

topology
new_topology (int nleaves) 
{
  topology tree;
  int i;
	size_t sizeof_node = sizeof (struct topol_node_struct);

  if (nleaves > MaxLeaves)
	biomc2_error ( "number of leaves > %d\n", MaxLeaves);

  tree = 
  (topology) biomc2_malloc (sizeof (struct topology_struct));
  tree->nleaves = nleaves;
  tree->nnodes  = 2*nleaves - 2;
  tree->root = NULL;
	tree->prev = tree->next = NULL;
	tree->first_segment = tree->last_segment = -1;
	tree->likelihood_accepted = tree->likelihood_current = tree->likelihood_proposal = 0.;

  tree->nodelist =
  (topol_node*) biomc2_malloc (tree->nnodes * sizeof (topol_node));

  /* tree->nodelist will store the actual nodes */
  for (i=0; i<tree->nleaves; i++) { 
    tree->nodelist[i] = (topol_node) biomc2_malloc (sizeof_node);
		tree->nodelist[i]->internal = false;
		tree->nodelist[i]->bipartition.node = tree->nodelist[i];
		split_initialize_l (&(tree->nodelist[i]->bipartition), i);
		tree->nodelist[i]->u_done = false; 
		tree->nodelist[i]->d_done = true;
		tree->nodelist[i]->left = tree->nodelist[i]->right = NULL;
		tree->nodelist[i]->map_id = tree->nodelist[i]->id = i;
  }
	tree->nodelist[0]->u_done = true;
	tree->nodelist[0]->d_done = false;
  for (i=tree->nleaves; i<tree->nnodes; i++) { 
    tree->nodelist[i] = (topol_node) biomc2_malloc (sizeof_node);
		tree->nodelist[i]->internal = true;
		tree->nodelist[i]->bipartition.node = tree->nodelist[i];
		tree->nodelist[i]->u_done = tree->nodelist[i]->d_done = true;
		tree->nodelist[i]->map_id = tree->nodelist[i]->id = i;
  }

	tree->root = tree->nodelist[0];
	tree->root->up = tree->root->right = NULL;

  return tree;
}

void 
del_topology (topology tree) 
{
  if (tree) {
    if (tree->nodelist) {
      int i;
      for (i=tree->nnodes-1; i >=0; i--) {
        if (tree->nodelist[i]) free (tree->nodelist[i]);
      }
      free (tree->nodelist);
    }
    free (tree);
  }
}

void
copy_topology_from_nexus_tree (topology tree, nexus_tree nxs_tree)
{
  int i, id, node_id;

  for (i = 0; i < tree->nnodes; i++) {
    node_id = nxs_tree->nodelist[i]->id;
		tree->nodelist[node_id]->map_id = tree->nodelist[node_id]->id = node_id;

    if (nxs_tree->nodelist[i]->up) {
      id = nxs_tree->nodelist[i]->up->id;
      tree->nodelist[node_id]->up = tree->nodelist[id];
    }
		else tree->nodelist[node_id]->up = NULL; 
		
		if (nxs_tree->nodelist[i]->left) {
			id = nxs_tree->nodelist[i]->left->id;
			tree->nodelist[node_id]->left = tree->nodelist[id];
		}
		else tree->nodelist[node_id]->left = NULL; 
		
		if (nxs_tree->nodelist[i]->right) {
			id = nxs_tree->nodelist[i]->right->id;
      tree->nodelist[node_id]->right = tree->nodelist[id];
		}
    else tree->nodelist[node_id]->right = NULL;
  } // for (nnodes)

  id = nxs_tree->root->id;
  if (id) biomc2_error ( "root id in nexus_tree is not zero\n");
  tree->root = tree->nodelist[id];
	
	for (i=0; i<tree->nleaves; i++) {
		tree->nodelist[i]->u_done = tree->nodelist[i]->d_done = true;
  }
	for (i=tree->nleaves; i<tree->nnodes; i++) {
		tree->nodelist[i]->u_done = tree->nodelist[i]->d_done = false;
	}

	update_topology_sisters (tree);
}


void
copy_topology_from_topology_mapping (topology to_tree, topology from_tree, split_space split)
{
	int i;
	for (i=from_tree->nleaves; i<from_tree->nnodes; i++) 
		from_tree->nodelist[i]->u_done = from_tree->nodelist[i]->d_done = false;

	/* from_tree receives u_node, d_done info (up and down subtrees are the same
	 * or not) */
	update_topology_from_topology_split (from_tree, to_tree, split);

	/* u_done, d_done info is copied to to_tree and original from_tree values are
	 * restored (they should be all TRUE since the topology was accepted in MCMC) */
	for (i=to_tree->nleaves; i < to_tree->nnodes; i++) {
		to_tree->nodelist[i]->u_done = from_tree->nodelist[i]->u_done;
		to_tree->nodelist[i]->d_done = from_tree->nodelist[i]->d_done;
		from_tree->nodelist[i]->u_done = from_tree->nodelist[i]->d_done = true;
	}

	copy_topology_from_topology (to_tree, from_tree);

}

void
copy_topology_from_topology (topology to_tree, topology from_tree)
{
	int i;

	to_tree->nodelist[0]->left = to_tree->nodelist[from_tree->nodelist[0]->left->id];
	for (i=1; i < from_tree->nleaves; i++) {
		to_tree->nodelist[i]->up = to_tree->nodelist[from_tree->nodelist[i]->up->id];
	}

  for (i = from_tree->nleaves; i < from_tree->nnodes; i++) {
//    if (from_tree->nodelist[i]->id != i) 
//    biomc2_error ("from %d to %d [%d not equal %d]\n", from_tree->id, to_tree->id, from_tree->nodelist[i]->id, i);
    to_tree->nodelist[i]->map_id = from_tree->nodelist[i]->map_id;

		to_tree->nodelist[i]->up    = to_tree->nodelist[from_tree->nodelist[i]->up->id];
		to_tree->nodelist[i]->left  = to_tree->nodelist[from_tree->nodelist[i]->left->id];
		to_tree->nodelist[i]->right = to_tree->nodelist[from_tree->nodelist[i]->right->id];
  } // for (nnodes)

  update_topology_sisters (to_tree);
}

void 
update_topology_sisters (topology tree)
{
	int i;

	/* skips root leaf */ 
	for (i=1; i < tree->nnodes; i++) {
    /* ingroup root is left child and has no sister */
		if (tree->nodelist[i]->up->left == tree->nodelist[i])
			tree->nodelist[i]->sister = tree->nodelist[i]->up->right;
		else
			tree->nodelist[i]->sister = tree->nodelist[i]->up->left;
	}
}

bool 
node1_is_child_of_node2 (topol_node node1, topol_node node2)
{
  if (node1 == node2) return true;
  if (node1->up && node1_is_child_of_node2 (node1->up, node2)) return true;
  return false;
}

topol_node
reroot_topology (topology tree, topol_node newroot) 
{
	if (newroot->internal)
		biomc2_error ( "cannot root internal node\n");
	if (newroot != tree->root) {
		topol_node pup = newroot->up, p = newroot, tmp = newroot->up->up;
		newroot->left = pup;
		newroot->up = NULL;
		while (pup != tree->root)
		 {                       /* p -> pup -> tmp (tmp may be old root) */
			if (pup->left == p)
				pup->left = tmp;
			else
				pup->right = tmp;
			pup->up = p;
			p = pup;
			pup = tmp;
			tmp = tmp->up;
		 }
		/* p -> pup (=oldroot) */
		pup->left = NULL;
		pup->up = p;
	}
	tree->root = newroot;
	return newroot;
}

char *
topology_to_string (const topology tree) 
{
	char *str;
	/* allocate space for str (overestimate size) */
	int size = tree->nnodes * 4 + tree->nleaves * 8;
	
	str = (char *) biomc2_malloc (sizeof (char) * size);
  memset (str, 0, sizeof (char) * size);
	sprintf (str, "(%d,", tree->root->id+1);
	topology_subtree_to_string (str, tree->root->left);
  sprintf (str, "%s);", str);
	return str;
}

void
topology_subtree_to_string (char *str, const topol_node node)
{
	if (node->internal) { /* internal node */
		sprintf (str, "%s(", str);
		topology_subtree_to_string (str, node->left);
		sprintf (str, "%s,", str);
		topology_subtree_to_string (str, node->right);
		sprintf (str, "%s)", str);
	}
	else sprintf (str, "%s%d", str, node->id+1);
}


char *
topology_to_string_with_name (const topology T, char **taxlabel)
{
 	char *str;
  int size = 1 + T->nnodes * 4, i;
  /* allocate space for str (overestimate size) */
  for (i=0; i < T->nleaves; i++) size += strlen (taxlabel[i]) + 1;

  str = (char *) biomc2_malloc (sizeof (char) * size);
  memset (str, 0, sizeof (char) * size);
  sprintf (str, "(%s,", taxlabel[T->root->id]);
	topology_subtree_to_string_with_name (str, T->root->left, taxlabel);
  sprintf (str, "%s);",str);
  return str;
}

void
topology_subtree_to_string_with_name (char *str, const topol_node node, char **taxlabel)
{
  if (node->internal) { /* internal node */
    sprintf (str, "%s(", str);
    topology_subtree_to_string_with_name (str, node->left, taxlabel);
    sprintf (str, "%s,", str);
		topology_subtree_to_string_with_name (str, node->right, taxlabel);
    sprintf (str, "%s)", str);
  }
  else sprintf (str, "%s%s", str, taxlabel[node->id]);
}

char *
topology_to_string_with_name_generate_length (const topology T, char **taxlabel, double min_length, double max_length)
{
  char *str;
  int size = 1 + T->nnodes * 4, i;
  /* allocate space for str (overestimate size) */
  for (i=0; i < T->nleaves; i++) size += strlen (taxlabel[i]) + 1;
  size *= 6;

  str = (char *) biomc2_malloc (sizeof (char) * size);
  memset (str, 0, sizeof (char) * size);
	sprintf (str, "(%s:%.4f,", taxlabel[T->root->id], 
					 min_length + max_length * biomc2_rng_uniform_pos (random_number));
  topology_subtree_to_string_with_name_generate_length (str, T->root->left, taxlabel, &min_length, &max_length);
  sprintf (str, "%s):%.4f;",str, min_length + max_length * biomc2_rng_uniform_pos (random_number));
  return str;
}

void
topology_subtree_to_string_with_name_generate_length (char *str, const topol_node node, char **taxlabel, 
																											double *min_length, double *max_length)
{
  if (node->internal) { /* internal node */
    sprintf (str, "%s(", str);
    topology_subtree_to_string_with_name_generate_length (str, node->left, taxlabel, min_length, max_length);
    sprintf (str, "%s,", str);
    topology_subtree_to_string_with_name_generate_length (str, node->right, taxlabel, min_length, max_length); 
    sprintf (str, "%s):%.4f", str, (*min_length) + (*max_length) * biomc2_rng_uniform_pos (random_number));
  }
	else sprintf (str, "%s%s:%.4f", str, taxlabel[node->id], 
								(*min_length) + (*max_length) * biomc2_rng_uniform_pos (random_number));
}


