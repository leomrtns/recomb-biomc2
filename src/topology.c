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


topology
new_topology (int nleaves) 
{
  topology tree;
  int i;
	size_t sizeof_node = sizeof (struct topol_node_struct);



  if (nleaves < 3) biomc2_error ("an unrooted tree needs at least three leaves");

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
		tree->nodelist[i]->split = new_bipartition (nleaves);
		bipartition_initialize (tree->nodelist[i]->split, i);
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
		tree->nodelist[i]->split = new_bipartition (nleaves);
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
        del_bipartition (tree->nodelist[i]->split);
        free (tree->nodelist[i]);
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

/* Return the next free character so emitting a tree takes linear time in its Newick length. */
static char *
write_subtree (char *out, topol_node node, char **labels, const double *lengths)
{
  if (node->internal) {
    *out++ = '(';
    out = write_subtree (out, node->left, labels, lengths);
    *out++ = ',';
    out = write_subtree (out, node->right, labels, lengths);
    *out++ = ')';
  }
  else if (labels) out += sprintf (out, "%s", labels[node->id]);
  else out += sprintf (out, "%d", node->id + 1);
  if (lengths) out += sprintf (out, ":%.4f", lengths[0] + lengths[1] * biomc2_rng_uniform_pos (random_number));
  return out;
}

static char *
write_tree (topology tree, char **labels, const double *lengths)
{
  size_t size = (size_t)tree->nnodes * 4 + 8;
  for (int i=0; i < tree->nleaves; i++) size += labels ? strlen (labels[i]) : 12;
  if (lengths) size += ((size_t)tree->nnodes + 1) * 320; /* decimal expansion of a finite double */
  char *result = biomc2_malloc (size), *out = result;
  *out++ = '(';
  out = write_subtree (out, tree->root, labels, lengths);
  *out++ = ',';
  out = write_subtree (out, tree->root->left, labels, lengths);
  *out++ = ')';
  if (lengths) out += sprintf (out, ":%.4f", lengths[0] + lengths[1] * biomc2_rng_uniform_pos (random_number));
  strcpy (out, ";");
  return result;
}

char *
topology_to_string (const topology tree)
{
  return write_tree (tree, NULL, NULL);
}

char *
topology_to_string_with_name (const topology tree, char **labels)
{
  return write_tree (tree, labels, NULL);
}

char *
topology_to_string_with_name_generate_length (const topology tree, char **labels, double min, double max)
{
  double lengths[] = {min, max};
  return write_tree (tree, labels, lengths);
}
