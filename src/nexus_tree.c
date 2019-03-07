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

/*! \file nexus_tree.c
 *  \ingroup nexus 
 *  \brief Reads tree files in nexus format.
 */

#include "nexus.h"

#define DEFAULTBLENGTH 0. /*!< \brief Default branch length. */

/* local funtion prototypes */

/*! \brief Allocates memory for nexus_tree_struct. */
nexus_tree new_nexus_tree (int nleaves);

/*! \brief Allocates memory for nexus_treespace_struct (set of trees present in
 * nexus file).  */
nexus_treespace new_nexus_treespace (void);

/*! \brief Frees memory used by tree. */
void del_nexus_tree (nexus_tree T);

/*! \brief Reads tree from file and store in treespace. */
void add_tree_to_nexus_treespace (nexus_treespace tsp, const char *string, bool translate);

/*! \brief Reads translation table (one line) of the form "number = taxa name" in tree file. */
void translate_taxa_nexus_treespace (nexus_treespace tsp, char *string);

/*! \brief Creates nexus_tree structure. 
 *  
 * Given a string with a tree in newick format, initialize the tree and reroot
 * it on an arbitrary leaf (the first appearing on the string).
 * \param[in] string string with tree information in newick format (with or
 * without branch lengths)
 * \return newly created nexus_tree_struct
 */
nexus_tree new_nexus_tree_from_string (char **string);

/*! \brief Recursive function that creates a node based on parenthetic structure. */
nexus_node subtree_nexus_tree (nexus_tree tree, char *lsptr, char *rsptr, int *node_id, nexus_node up);

/*! \brief Reads leaf name (or number, if translation table is present). */
char* read_taxlabel ( const char *name_start, const char *name_end);

/*! \brief Preorder initialization of leaves. */
void create_leaflist_nexus_tree (nexus_tree tree, nexus_node this, int *id);

/*! \brief Remove the original root node (useless) and reroot on an arbitrary leaf. */
nexus_node unroot_nexus_tree (nexus_node oldroot, nexus_node newroot);

/*! \brief Counts the number of leaves and resolves (one) trifurcation of tree string. */
int number_of_leaves_in_newick (char **string);

/*! \brief Preorder initialization of _internal_ nodes; 'id' should be >= nleaves. */
void create_node_id_nexus_tree (nexus_node this, int *id);

/*! \brief Searches for (last reported) branch length on string or return default value.  */
double read_branch_length (char *right_string_ptr);

/*! \brief Returns position of innermost comma (divides string into two subtrees). */
int find_branch_split_newick (char *left_string_ptr, char *right_string_ptr);

/*! \brief Recursive function that prints node info to string in top-bottom order. */
void nexus_subtree_to_string (char *str, const nexus_node node, const bool has_branches);

nexus_tree
new_nexus_tree (int nleaves) 
{
	nexus_tree tree;
	int i, sizeof_node = sizeof (struct nexus_node_struct);

	tree = 
	(nexus_tree) biomc2_malloc (sizeof (struct nexus_tree_struct));
	tree->root = (nexus_node) biomc2_malloc (sizeof (struct nexus_node_struct));
	tree->nleaves = nleaves;
	tree->nnodes  = 2*nleaves - 2;
 	tree->has_branches = false;
	
	tree->nodelist =
	(nexus_node*) biomc2_malloc (tree->nnodes * sizeof (nexus_node));
	tree->leaflist =
	(nexus_node*) biomc2_malloc (tree->nleaves * sizeof (nexus_node));
	
  /* tree->nodelist will store the actual nodes */
  for (i=0; i<tree->nnodes; i++) { 
    tree->nodelist[i] = (nexus_node) biomc2_malloc (sizeof_node);
    tree->nodelist[i]->taxlabel = NULL;
  }

  return tree;
}

nexus_treespace
read_nexus_treespace_file (char *seqfilename)
{
  nexus_treespace treespace=NULL;
  FILE *seqfile;
  char *line=NULL, *needle_tip, *translate_taxa=NULL;
	bool option_begin_trees    = false,
			 option_translate_perm = false,
			 option_translate_temp = false; 
  size_t linelength = 0;


  seqfile = biomc2_fopen (seqfilename, "r");
  biomc2_getline (&line, &linelength, seqfile);

  if (!strcasestr (line, "NEXUS")) biomc2_error ( "Not a Nexus tree file\n");

	/* the variable *line should point always to the same value (no line++ or
	 * alike) */
	while (biomc2_getline (&line, &linelength, seqfile) != -1) {
		line = remove_nexus_comments (&line, &linelength, seqfile);

		/* do not parse empty lines ( = "\n\0") */	
		if (strlen (line) > 1) { 
			/* don't do anything untill 'BEGIN TREES' block */
      if ((!option_begin_trees) && (strcasestr (line, "BEGIN TREES"))) {
        option_begin_trees = true;
        treespace = new_nexus_treespace();
      } 

			/* check if we need to translate; in any case see if we have trees to read */
      else if (!option_translate_temp) {
        if (strcasestr (line, "TRANSLATE")) {
          option_translate_perm = true;
          option_translate_temp = true;
        }
        else if (strcasestr (line, "TREE") && 
								 (needle_tip = strcasestr (line, "="))) {
					needle_tip++; /* remove "=" from string */
					add_tree_to_nexus_treespace (treespace, needle_tip, option_translate_perm);
				}
			}
		
      /* we are reading translation table token <-> taxlabel */  
      if (option_translate_temp) {
        translate_taxa_nexus_treespace (treespace, line);
        if (strchr (line, ';')) option_translate_temp = false;
      }
      
		} // if (line)
	} //while (biomc2_getline)
	
	fclose (seqfile);
	if (translate_taxa) free (translate_taxa);
	if (line) free (line);
	
	return treespace;
}

nexus_treespace
new_nexus_treespace (void) 
{
	nexus_treespace tsp;
	
	tsp = (nexus_treespace) biomc2_malloc (sizeof (struct nexus_treespace_struct));

	tsp->ntrees  = 0;
	tsp->nleaves = 0;
	tsp->nnodes  = 0;
  tsp->T = NULL;
  tsp->taxlabel = NULL;
	/* tsp->T vector is increased by add_tree_to_nexus_treespace() 
	 * tsp->taxlabel vector is setup by translate_taxa_nexus_treespace() or
	 * add_tree_to_nexus_treespace(), whichever is used
	 * tsp->taxlabel_hash is setup by one of the above */
	return tsp;
}

void
del_nexus_treespace (nexus_treespace tsp) 
{
  if (tsp) {
    int i;
		if (tsp->T) {
      for (i=tsp->ntrees-1; i>=0; i--) if (tsp->T[i]) {
        del_nexus_tree (tsp->T[i]);
			}
			free (tsp->T);
		}
		if (tsp->taxlabel) { 
			for (i=tsp->nleaves-1; i>=0; i--) 
				if (tsp->taxlabel[i])  free (tsp->taxlabel[i]);
			free (tsp->taxlabel);
		}
		del_hashtable (tsp->taxlabel_hash);
		free (tsp);
  }
}

void 
del_nexus_tree (nexus_tree T) 
{
	if (T) {
		int i;
		for (i=T->nnodes-1; i >=0; i--) {
			if (T->nodelist[i]) {
				/* taxlabels cannot be freed here since they are shared in treespace 
				 * COMMENT from BUG 2008.08.14: if they are shared (pointers to singly initialized tsp->taxlabel) there is no
				 * need to free them */	
				/*	if (T->nodelist[i]->taxlabel) free (T->nodelist[i]->taxlabel); */
				free (T->nodelist[i]);
			}
		}
		if (T->nodelist) free (T->nodelist);
		if (T->leaflist) free (T->leaflist);
		free (T);
	}
}

void
add_tree_to_nexus_treespace (nexus_treespace tsp, const char *string, bool translate) 
{
  int i, index;
  char *local_string;
  nexus_tree tree;

  /* use local copy of string to avoid problems with biomc2_getline() */
  local_string = (char*) biomc2_malloc (sizeof (char) * (strlen (string) + 1));
  strcpy (local_string, string);
  
  tree	= new_nexus_tree_from_string (&local_string);
  if (local_string) free (local_string);

 /* first tree read and no TRANSLATE command in nexus file */
  if ((tsp->ntrees == 0) && (!translate)) {
    tsp->nnodes  = tree->nnodes;
    tsp->nleaves = tree->nleaves;
    
    /* tsp->taxlabel will point to names of first tree. This info will be also 
     * available at the hashtable */ 
    tsp->taxlabel = (char**) biomc2_malloc (sizeof (char*) * tsp->nleaves);
    tsp->taxlabel_hash = new_hashtable (tsp->nleaves);
    
    for (i=0; i< tsp->nleaves; i++) { 
      tsp->taxlabel[i] = tree->leaflist[i]->taxlabel;
      insert_hashtable (tsp->taxlabel_hash, tsp->taxlabel[i], i);
	  tree->leaflist[i]->id = i;
    }
	create_node_id_nexus_tree (tree->root->left, &i);
  }
  
  /* first tree read and TRANSLATE command in nexus file */
  else if ((tsp->ntrees == 0) && (translate)) {
    for (i=0; i< tsp->nleaves; i++) {
      sscanf (tree->leaflist[i]->taxlabel, " %d ", &index);
      index--;
      if (index<0 || index >=tsp->nleaves)
        biomc2_error ( "leaf name out of range (1...NTAX) \n");
      if (tree->leaflist[i]->taxlabel) free (tree->leaflist[i]->taxlabel); 
      tree->leaflist[i]->taxlabel = tsp->taxlabel[index];
	  tree->leaflist[i]->id = index;
    }
	create_node_id_nexus_tree (tree->root->left, &i);
  }
  
  /* not the first tree read and no TRANSLATE command in nexus file */
  else if ((tsp->ntrees > 0) && (!translate)) {
    if (tsp->nleaves != tree->nleaves) {
      del_nexus_tree (tree);
      del_nexus_treespace (tsp);
      biomc2_error ( "number of leaves disagrees between trees\n");
    }
    /* use hashtable to check if names are consistent and point all leaves to
     * taxlabel vector */
    for (i=0; i< tree->nleaves; i++) {
      index = lookup_hashtable (tsp->taxlabel_hash, tree->leaflist[i]->taxlabel);
	  if (index < 0) biomc2_error ( "leaf names disagree between trees\n");
      if (tree->leaflist[i]->taxlabel) free (tree->leaflist[i]->taxlabel); 
      tree->leaflist[i]->taxlabel = tsp->taxlabel[index];
	  tree->leaflist[i]->id = index;
		}
		create_node_id_nexus_tree (tree->root->left, &i);
	}
  
  /* not the first tree read and TRANSLATE command in nexus file */
  else {
    if (tsp->nleaves != tree->nleaves) {
      del_nexus_tree (tree);
      del_nexus_treespace (tsp);
      biomc2_error ( "number of leaves disagrees between tree and TRANSLATE command\n");
    }
    for (i=0; i< tree->nleaves; i++) {
      sscanf (tree->leaflist[i]->taxlabel, " %d ", &index);
      index--;
      if (index<0 || index >=tsp->nleaves)
        biomc2_error ( "leaf name out of range (1...NTAX) \n");
      if (tree->leaflist[i]->taxlabel) free (tree->leaflist[i]->taxlabel); 
      tree->leaflist[i]->taxlabel = tsp->taxlabel[index];
	  tree->leaflist[i]->id = index;
		}
		create_node_id_nexus_tree (tree->root->left, &i);
	}
  
  tsp->ntrees++;
	tsp->T =
	(nexus_tree*) biomc2_realloc ((nexus_tree*) tsp->T, sizeof (nexus_tree) * (tsp->ntrees));
	tsp->T[tsp->ntrees-1] = tree;
}

void
translate_taxa_nexus_treespace (nexus_treespace tsp, char *string) 
{
  int i, index;
  char *c, *s, *last, *comma_position, seqname[MAX_NAME_LENGTH]="";


  /* the file may have the first token<->taxon_name at the same line as
   * the "TRANSLATE" command. */
  if ( (s = strcasestr (string, "TRANSLATE")) ) s += strlen ("TRANSLATE");
  else s = (string);
  last = s + strlen (s) +1;

  /* one or more "token1 taxon_id1, token2 taxon_id2," */
  while ( (comma_position = strchr (s, ',')) && 
          (sscanf (s, " %d %s ,", &index, seqname) == 2) && 
          (s<last)) 
   {
    tsp->nleaves++;
    tsp->taxlabel =
    (char**) biomc2_realloc ((char**) tsp->taxlabel, sizeof (char*) * tsp->nleaves);
    index--; /* in nexus we have index = 1...NTAX */
    while ( (c = strrchr (seqname, ',')) ) *c = '\0'; /* remove trash at the end */
    while ( (c = strrchr (seqname, ';')) ) *c = '\0';
    i = strlen (seqname) + 1;

    if ((i<0) || (index<0))
      biomc2_error ( "unexpected leaf name/location in TRANSLATE command\n"); 
    
    tsp->taxlabel[index] = (char*) biomc2_malloc (sizeof (char) * i);
    strcpy (tsp->taxlabel[index], seqname);
    s = comma_position+1;
   }
  
  /* maybe only one "token taxon_id" (the last line, e.g.) */
  if (sscanf (s, " %d %s ", &index, seqname) == 2) {
    tsp->nleaves++;
    tsp->taxlabel =
    (char**) biomc2_realloc ((char**) tsp->taxlabel, sizeof (char*) * tsp->nleaves);
    index--; /* in nexus we have index = 1...NTAX */
    while ( (c = strrchr (seqname, ',')) ) *c = '\0'; /* remove trash at the end */
    while ( (c = strrchr (seqname, ';')) ) *c = '\0';
    i = strlen (seqname) + 1;
    
    if ((i<0) || (index<0))
      biomc2_error ( "unexpected leaf name/location in TRANSLATE command\n"); 
    
    tsp->taxlabel[index] = (char*) biomc2_malloc (sizeof (char) * i);
    strcpy (tsp->taxlabel[index], seqname);
  }
  
  /* when we finished reading the leaf names */ 
  if (strchr (string, ';')) {
    tsp->taxlabel_hash = new_hashtable (tsp->nleaves);
    for (i=0; i<tsp->nleaves; i++) 
      insert_hashtable (tsp->taxlabel_hash, tsp->taxlabel[i], i);
  }
}

double
average_rate_nexus_tree (nexus_tree T)
{
  int i, j = 0;
  double rate = 0.;

  for (i=0; i< T->nnodes; i++)
    if (T->nodelist[i]->branch_length>0) {
      rate += T->nodelist[i]->branch_length;
      j++;
    }
	if (rate>0)
		return rate/(double)(j);
	else return -1.;
}

void map_hashtable_order (int *order, hashtable hash, char **taxlabel, int size)
{
  int i;

  for (i=0; i < size; i++) {
	order[i] = lookup_hashtable (hash, taxlabel[i]);
	if (order[i] < 0) biomc2_error ( "tree label %s not found in sequence data\n", taxlabel[i]); 
  }
}

void order_nexus_treespace_id (nexus_treespace t, int *order)
{
  int i, j, root_id;

  for (i=0; i < t->ntrees; i++) {
	root_id = -1;
	for (j=0; j < t->nleaves; j++) {
	  t->T[i]->leaflist[j]->id = order[ t->T[i]->leaflist[j]->id ];
	  if (!t->T[i]->leaflist[j]->id) root_id = j;
	}
	if (root_id < 0) biomc2_error ( "problem in alignment X tree mapping\n");
	reroot_nexus_tree (t->T[i], t->T[i]->leaflist[root_id]);
  }
}

nexus_tree
new_nexus_tree_from_string (char **string) 
{
  char *lptr, *rptr;
  int id, nleaves;
	nexus_tree T; 
	
	*string = remove_space_from_string (*string);
	nleaves = number_of_leaves_in_newick (string); /* resolves trifurcation */
  T = new_nexus_tree (nleaves);
  if (strchr (*string, ':')) T->has_branches = true;
	
	/* begin & end of string */
	lptr = *string;	
	rptr = *string + strlen (*string) - 1;
	
	id = -1; /* This function does the actual creation of the tree */
	T->root = subtree_nexus_tree (T, lptr, rptr, &id, NULL);
	
	id = 0; /* vector of pointers to the tree leaves */
	create_leaflist_nexus_tree (T, T->root, &id); 

	/* reroot the tree on a leaf and remove extra (innermost) node */
	T->root = unroot_nexus_tree (T->root, T->leaflist[0]);
	
  return T;
}

nexus_node
subtree_nexus_tree (nexus_tree tree, char *lsptr, char *rsptr, int *node_id, nexus_node up) 
{
	nexus_node thisnode;
	
	if (*node_id<0) thisnode = tree->root;
	else thisnode = tree->nodelist[*node_id];	
	thisnode->up = up;
	thisnode->id = -1; 
	thisnode->branch_length = tree->has_branches ? read_branch_length (rsptr) : DEFAULTBLENGTH;
	thisnode->left = NULL;
	thisnode->right = NULL;
	thisnode->taxlabel = NULL;

	(*node_id)++;

	if (*lsptr == '(') { /* internal node */
		char *newend = rsptr;
		int comma_pos = find_branch_split_newick (lsptr, rsptr);
		
		thisnode->left = subtree_nexus_tree (tree, lsptr+1, lsptr+comma_pos-1, node_id, thisnode);
		
		while (newend != lsptr && newend != NULL && *newend != ')') newend--;
		if (newend == lsptr) newend = rsptr;
		
		thisnode->right = subtree_nexus_tree (tree, lsptr+comma_pos+1, newend-1, node_id, thisnode);
 	}

  else { /* leaf */
		thisnode->taxlabel = read_taxlabel (lsptr, rsptr);
  }
	
  return thisnode;
}

char *
read_taxlabel ( const char *name_start, const char *name_end) 
{
	size_t seqsize;
  size_t i;
	char *tmp, *label=NULL;
  tmp = (char*) name_start;
  while ((tmp <= name_end) && (*tmp != ',') && (*tmp != ')') && (*tmp != ':')) tmp++;
  seqsize = tmp - name_start;
  i = sizeof (char)*(seqsize+1);
  label = (char*) biomc2_malloc (i);
  label[0] = '\0';
//  sprintf (label, '\0');
	strncat (label, name_start, seqsize);
	return label;
}

void
create_leaflist_nexus_tree (nexus_tree tree, nexus_node this, int *id) 
{
	if (this->taxlabel != NULL)
		tree->leaflist[(*id)++] = this;
	else {
		if (this->left) create_leaflist_nexus_tree (tree, this->left, id);
		if (this->right) create_leaflist_nexus_tree (tree, this->right, id);
	}
}

void
create_node_id_nexus_tree (nexus_node this, int *id) 
{
  this->id = (*id)++;
  if ( (this->left)  && (this->left->id < 0)  ) create_node_id_nexus_tree (this->left, id);
  if ( (this->right) && (this->right->id < 0) ) create_node_id_nexus_tree (this->right, id);
}

double
read_branch_length (char *right_string_ptr) 
{
 	char *backwards = right_string_ptr;
 	double branch=0.;
	
	if ((*backwards == ')') || (*backwards == ',')) return DEFAULTBLENGTH;
	while (*backwards != ':') backwards--;
	sscanf (backwards, ": %lf", &branch);
	if (branch < 0.) return 0.;
	else return branch;
}

int
find_branch_split_newick (char *left_string_ptr, char *right_string_ptr) 
{
	int nLevel = 0;
	int comma_position = 0;
	char *treeptr = left_string_ptr;
	
	while ((*treeptr != ',') || (nLevel != 1)) {
		if (*treeptr == '(') nLevel++;
		if (*treeptr == ')') nLevel--;
		treeptr++;
		comma_position++;
    if (treeptr == right_string_ptr) return comma_position; // added 20160809, to avoid complain of unused *right var
	}
	return comma_position;
}

nexus_node
unroot_nexus_tree (nexus_node oldroot, nexus_node newroot) 
{
	if (!newroot->taxlabel)
		biomc2_error ( "cannot root internal node\n");
	else {
		nexus_node pup = newroot->up, p = newroot, tmp;
		double tmpd1 = p->branch_length, tmpd2 = 0.0;
		p->branch_length = 0.0;
		if (pup != oldroot) { /* climbing up the tree */
			newroot->left = pup;  /* we define that root has only left child */
			newroot->up = NULL;
			while (pup->up != oldroot) {
				/* (newroot) ...  p -> pup -> tmp  ... (oldroot) */
				tmp = pup->up;    /*otherwise we loose track of old parent */
				if (pup->left == p)
					pup->left = tmp;
				else
					pup->right = tmp;
				pup->up = p;
				tmpd2 = pup->branch_length;
				pup->branch_length = tmpd1;
				tmpd1 = tmpd2;
				p = pup;
				pup = tmp;
			}                   // while ()
			/* now pup->up points to old root, that will be removed */
			/* p -> pup -> oldroot <- tmp */
			if (pup->up->left == pup)
				tmp = pup->up->right;
			else
				tmp = pup->up->left;
			if (pup->left == p)
				pup->left = tmp;
			else
				pup->right = tmp;
			pup->up = p;
			tmp->up = pup;
			tmp->branch_length += pup->branch_length;
			pup->branch_length = tmpd1;
		} //if (pup != oldroot)
		
		else { /* pup (=newroot->up) == oldroot */
			if (pup->left == p)
				newroot->left = pup->right;
			else
				newroot->left = pup->left;
			newroot->up = NULL;
			newroot->left->up = newroot;
			newroot->left->branch_length += tmpd1;
		}
		free (oldroot);
	} // else() (not internal node)
  return newroot;
}

nexus_node
reroot_nexus_tree (nexus_tree tree, nexus_node newroot) 
{
 	if (!newroot->taxlabel)
 		biomc2_error ( "cannot root internal node\n");
	if (newroot != tree->root) {
		nexus_node pup = newroot->up, p = newroot, tmp = newroot->up->up;
		double tmpd1 = newroot->branch_length, tmpd2 = 0.0;
		newroot->branch_length = 0.0;
		newroot->left = pup;
		newroot->up = newroot->right = NULL;
		while (pup != tree->root)
 		 {                       /* p -> pup -> tmp (tmp may be old root) */
			if (pup->left == p)
				pup->left = tmp;
			else
				pup->right = tmp;
			pup->up = p;
			tmpd2 = pup->branch_length;
			pup->branch_length = tmpd1;
			tmpd1 = tmpd2;
			p = pup;
			pup = tmp;
			tmp = tmp->up;
		 }
		/* p -> pup (=oldroot) */
		pup->left = NULL;
		pup->up = p;
		pup->branch_length = tmpd1;
	}
	tree->root = newroot;
	return newroot;
}

int
number_of_leaves_in_newick (char **string) 
{
	int nopen = 0, nclose = 0, ncommas = 0, i, nsplit = 0;
	int len = strlen (*string);
	bool has_branches = false;
	if (*(*string + len - 1) == ';')
		*(*string + len - 1) = '\0';
	for (i = 0; i < len; i++) {
		char current = (*string)[i];
		if (current == ',' && (nopen - nclose) == 1) {
			if (!nsplit) nsplit = i;
			ncommas++;
		}
		else if (current == '(')
			nopen++;
		else if (current == ')')
			nclose++;
		else if (current == ':')
			has_branches = true;
	}
	
	if (nopen != nclose || ncommas > 2 || ncommas < 1)
		biomc2_error ( "Invalid tree structure: %s", *string);
	/* try to root an unrooted tree */
	if (ncommas == 2)
 	 {
		char *tstring = NULL;
		int add = 3;
		if (has_branches) add = 7;
		*string =
		(char *) biomc2_realloc ((char *) *string, sizeof (char) * (len + add + 1));
		tstring = 
		(char *) biomc2_malloc (sizeof (char) * (len + add + 1));
		bzero (tstring, sizeof (char) * (len + add + 1));
		tstring = strncat (tstring, *string, nsplit);
		tstring = strcat (tstring, ",(");
		tstring = strncat (tstring, *string + nsplit + 1, len - nsplit - 2);
		if (has_branches)
			tstring = strcat (tstring, "):0.0");
		else
			tstring = strcat (tstring, ")");
		strcpy (*string, tstring);
		free (tstring);
		nopen++;
	 }
	return nopen + 1;
}

char *
nexus_tree_to_string (const nexus_tree tree) 
{
	int size = 1, i;
	char *str;
	/* allocate space for str (overestimate size) */
	for (i = 0; i < tree->nleaves; i++)
		size += strlen (tree->leaflist[i]->taxlabel) + 1;
	size += tree->nnodes * 4;
	if (tree->has_branches)
		size += tree->nnodes * 12;
	
	str = (char *) biomc2_malloc (sizeof (char) * size);
  memset (str, 0, sizeof (char) * size);
  sprintf (str, "(");
	nexus_subtree_to_string (str, tree->root->left, tree->has_branches);
	if (tree->has_branches)
		sprintf (str, "%s,%s: %.8f);", str, tree->root->taxlabel, tree->root->branch_length);
	else
		sprintf (str, "%s,%s);", str, tree->root->taxlabel);
	return str;
}

void
nexus_subtree_to_string (char *str, const nexus_node node, 
															const bool has_branches)
{
	if (node->taxlabel != NULL) {
		if (has_branches)
			sprintf (str, "%s%s: %.8f", str, node->taxlabel, node->branch_length);
		else
			sprintf (str, "%s%s", str, node->taxlabel);
	}
	else
	 { /* internal node */
	 	sprintf (str, "%s(", str);
		nexus_subtree_to_string (str, node->left, has_branches);
		sprintf (str, "%s,", str);
		nexus_subtree_to_string (str, node->right, has_branches);
		if (has_branches)
			sprintf (str, "%s): %.8f", str, node->branch_length);
		else
			sprintf (str, "%s)", str);
	 }
}

void
graphviz_file_nexus_tree (FILE * fout, char *label, const nexus_tree tree) 
{
 	int i;
 	fprintf (fout, "graph G {\n");
 	fprintf (fout, "  graph [ size=\"7,9\" page=\"8.5,11\" center=\"\" ]\n");
 	fprintf (fout, "  node  [ fontsize = \"8\" width=.08, hight=.08 ]\n");
 	fprintf (fout, "  edge  [ fontsize = \"6\" len=1.5 ]\n");
	for (i = 0; i < tree->nnodes; i++)
	 {
		if (tree->nodelist[i]->taxlabel)
			fprintf (fout,
							 "  %d\t[ label = \"%d\\n%s\" width=.16, hight=.16 ];\n",
							 tree->nodelist[i]->id, tree->nodelist[i]->id,
							 tree->nodelist[i]->taxlabel);
		if (tree->nodelist[i]->left)
			fprintf (fout, "  %d -- %d\t[ label = \"%f\" ];\n",
							 tree->nodelist[i]->id, tree->nodelist[i]->left->id,
							 tree->nodelist[i]->left->branch_length);
		if (tree->nodelist[i]->right)
			fprintf (fout, "  %d -- %d\t[ label = \"%f\" ];\n",
							 tree->nodelist[i]->id,
							 tree->nodelist[i]->right->id,
							 tree->nodelist[i]->right->branch_length);
	 }
	if (label)
		fprintf (fout, "  label =\"%s\";\n", label);
	fprintf (fout, "\n}\n");
	fflush (fout);
}

