/*****************************************************************************
#   Copyright (C) 1993-2009 by Phil Green.                          
#   All rights reserved.                           
#                                                                           
#   This software is part of a beta-test version of the swat/cross_match/phrap 
#   package.  It should not be redistributed or
#   used for any commercial purpose, including commercially funded
#   sequencing, without written permission from the author and the
#   University of Washington.
#   
#   This software is provided ``AS IS'' and any express or implied
#   warranties, including, but not limited to, the implied warranties of
#   merchantability and fitness for a particular purpose, are disclaimed.
#   In no event shall the authors or the University of Washington be
#   liable for any direct, indirect, incidental, special, exemplary, or
#   consequential damages (including, but not limited to, procurement of
#   substitute goods or services; loss of use, data, or profits; or
#   business interruption) however caused and on any theory of liability,
#   whether in contract, strict liability, or tort (including negligence
#   or otherwise) arising in any way out of the use of this software, even
#   if advised of the possibility of such damage.
#                                                                         
#*****************************************************************************/


/* routines for dealing with sequence files in FASTA format */

#include "swat.h"

#define ID_BUFFERSIZE 50 /* max size of ID field in database entries -- is adjusted upward at runtime if nec .*/
#define DESCRIP_BUFFERSIZE 100
#define MAX_NUM_BLOCKS 5000 
#define BLOCK_SIZE 10000 /* was 1000 -- increased 10/5/07 */

extern Parameters *parameters;

unsigned char area_comp_mat[256];
int qseq_trans[256];
int t_num_entries;

typedef long int POS;
static POS id_pos, descrip_pos, seq_pos, seq_length;
/* N.B. following needed by Bonfield code, so original "static" designation has
   been removed */
Database *head_db;
Seq_entry *seq_block_table[MAX_NUM_BLOCKS];
int num_blocks;
unsigned char area_alph_mat[256], buffer_comp_mat[256], buffer_alph_mat[256]; /* one version for seq_area, another for seq_buffer */

/* called in parameters.c. 

 sets area_alph_mat (used in alloc_mem and file_to_mem): 0 except for alphabet chars or replacements, preserves those or maps
       to upper case. Not needed for calf files; with compact_qual (internal conversion to calf) sets to base/qual byt
       with non-N, non-ACGT (or lower case), non-'.' going to 128

 buffer_alph_mat (used in  get_next_file_entry, get_seq, get_comp_seq)  0 except for alphabet chars or replacements, preserves those or maps
       to upper case.

 area_comp_mat (used in store_complements, and in append_diffdata in qual.c)  complements characters in the seq area.

 buffer_comp_mat (used in get_comp_seq) 

*/

set_mats(allow_buffer_lc, allow_area_lc)
     int allow_buffer_lc, allow_area_lc;
{
  int i, c, d, q;
  char vec[20], comp_vec[20];

  for (i = 0; i < 256; i++) {
    qseq_trans[i] = i;
    area_alph_mat[i] = !isalpha(i) ? 0 : (allow_area_lc ? i : toupper(i));
    buffer_alph_mat[i] = !isalpha(i) ? 0 : (allow_buffer_lc ? i : toupper(i));
  }

  if (parameters->DNA_flag) {
    area_alph_mat['.'] = 'N';
    buffer_alph_mat['-'] = buffer_alph_mat['.'] = 'N';

    strcpy(vec, "ACGTBDHVMRWSYK");
    strcpy(comp_vec, "TGCAVHDBKYWSRM");

    for (i = 0; i < 256; i++) buffer_comp_mat[i] = i;

    for (i = 0; c = vec[i]; i++) {
      buffer_comp_mat[c] = comp_vec[i];
      buffer_comp_mat[tolower(c)] = tolower(comp_vec[i]);
    }

    if (parameters->DNA_flag >= 3) {
      for (i = 0; i < 256; i++) {
	q = i & 63;
	area_comp_mat[i] = !q || q > 61 ? i : (192 - (i & 192)) + q;
	qseq_trans[i] = q ? "ACGT"[i >> 6] : "!N-*"[i >> 6];
      }
      if (parameters->DNA_flag == 3) {
	area_alph_mat['-'] = '-'; /* will be converted to 128 */
	d = parameters->default_qual;
	if (d > 60 || d < 0) fatalError("default qual must be <= 60 and >= 0");
	d++; /* add one, to accord with CALF format */
	for (i = 0; i < 256; i++) {
	  if (c = area_alph_mat[i])
	    area_alph_mat[i] = c == 'A' ? d : c == 'C' ? d + 64 : c == 'G' ? d + 128 : c == 'T' ? d + 192 
	      : c == 'N' ? 64 : 128; /* preserves qual 0; also need reserved char for X / stop */
	  /* note assumes letters are upper-case */
	}
      }
    }
    else {
      area_alph_mat['-'] = 'N'; /* CHECK THIS */
      for (i = 0; i < 256; i++) area_comp_mat[i] = i;
      for (i = 0; c = vec[i]; i++) {
	area_comp_mat[c] = i < 8 || !allow_area_lc ? comp_vec[i] : tolower(comp_vec[i]);
	area_comp_mat[tolower(c)] = i < 8 || !allow_area_lc  ? tolower(comp_vec[i]) : comp_vec[i];
     }
    }
  }
}

Database *append_db(file, in_memory, store_comps, store_aligns)
     File *file;
     int in_memory, store_comps, store_aligns;
{
  char *our_alloc();
  Database *db;

  db = (Database *)our_alloc(sizeof(Database));
  db->file = file;
  file->db = db;
  db->complements = store_comps;
  db->in_memory = in_memory;

  db->t_length = db->num_entries = 0;
  db->first_entry = -1;
  db->last_entry = -2;
  db->id_area = db->descrip_area = 0;
  db->orig_qual_area = db->adj_qual_area = 0;
  db->seq_area = 0;
  db->id_area_size = db->descrip_area_size = db->seq_area_size = 0;
  db->align_array = 0;

  db->id_buffer_size = ID_BUFFERSIZE;
  db->descrip_buffer_size = DESCRIP_BUFFERSIZE;
  
  db->seq_buffer_size = db->comp_seq_buffer_size = 0; /* SEQ_BUFFERSIZE; */
  db->seq_buffer = db->comp_seq_buffer = 0;

  db->id_buffer = (char *)our_alloc((db->id_buffer_size + 1) * sizeof(char));
  db->descrip_buffer = (char *)our_alloc((db->descrip_buffer_size + 1) * sizeof(char));
  /*
  db->seq_buffer = (char *)our_alloc((db->seq_buffer_size + 1) * sizeof(char));
  db->comp_seq_buffer = (char *)our_alloc((db->comp_seq_buffer_size + 1) * sizeof(char));
  db->comp_seq_buffer[0] = 0;
  */

  db->id_buffer[0] = db->descrip_buffer[0] = 0;
  db->seq_buffer_entry = db->comp_seq_buffer_entry = -1;

  db->next = head_db;
  head_db = db;

  if (in_memory) {
    notify("Reading query file into memory ...");
    alloc_mem(db, store_comps);
    file_to_mem(db);
    notify(" Done\n");

    if (store_comps && parameters->DNA_flag != 3) { /* DNA_flag = 3 case done after reading quals instead */
      store_complements(db); /* should really do later instead, after reading quals, for all DNA_flag vals -- but might break gcphrap or gccross_match */
    }
    if (store_aligns) {
      notify("Allocating align_entries ...");

      init_aligns(db);

      init_qual_areas(db);

      notify(" Done\n");
    }
  }
  return db;
}

store_complements(db)
     Database *db;
{
  SEQ_AREA j, k;
  unsigned char *seq_area;

  notify("Complementing ...");
  seq_area = db->seq_area;
  for (j = 0, k = db->seq_area_size - 1; j < k; j++, k--) 
    seq_area[k] = area_comp_mat[seq_area[j]];
  notify(" Done\n");

}

free_seq_buffers(db)
     Database *db;
{
  if (db->seq_buffer) {
    our_free(db->seq_buffer - 1);
    db->seq_buffer = 0;
    db->seq_buffer_size = 0;
  }
  if (db->comp_seq_buffer) {
    our_free(db->comp_seq_buffer - 1);
    db->comp_seq_buffer = 0;
    db->comp_seq_buffer_size = 0;
  }
  db->seq_buffer_entry = db->comp_seq_buffer_entry = -1;
}

init_aligns(db)
     Database *db;
{
  char *our_alloc();
  int i;
  Align_info *align_array;

  align_array = db->align_array 
    = (Align_info *)our_alloc(db->num_entries * sizeof(Align_info)) - db->first_entry;
  for (i = db->first_entry; i <= db->last_entry; i++) {

    align_array[i].seq_entry = align_array[i].equiv_class = i;

    align_array[i].template_start = align_array[i].template_end = 0;

    align_array[i].orig_qual = align_array[i].adj_qual = 0;

    align_array[i].segments = 0;

    align_array[i].first_start = align_array[i].rev_first_start = align_array[i].first_vec 
      = align_array[i].qual_start = get_seq_length(i); 
    align_array[i].last_end = align_array[i].rev_last_end = 0; /* should be -1 ?? */
    align_array[i].last_vec = align_array[i].qual_end = -1;

    align_array[i].anomalies = 0;
    align_array[i].blocked = 0;
    align_array[i].bypassed = 0;
    align_array[i].chemistry = 0;
    align_array[i].direction = 0;
    align_array[i].chimera_bits = 0;
    align_array[i].tags = 0;
  }
}

init_qual_areas(db)
     Database *db;
{
  char *our_alloc();
  Align_info *align_entry;
  SEQ_AREA i, t_size;
  int i_entry, incr;

  t_size = db->t_length;
  if (parameters->DNA_flag >= 3) {
    db->orig_qual_area = (char *)(db->seq_area + 1); 
    incr = 1;
    t_size += db->num_entries;
  }
  else {
    db->orig_qual_area = (char *)our_alloc(db->t_length * sizeof(char));
    for (i = 0; i < db->t_length; i++) db->orig_qual_area[i] = parameters->default_qual;
    incr = 0;
  }

/* with cross_match, adj_qual should be identical with orig_qual; otherwise, 
   make it separate */
  if (strcmp(parameters->calling_program, "cross_match")) {
    db->adj_qual_area = (char *)our_alloc(db->t_length * sizeof(char));
    for (i = 0; i < db->t_length; i++) db->adj_qual_area[i] = 0;
  }
  else db->adj_qual_area = db->orig_qual_area;

  for (i_entry = db->first_entry, i = 0; i_entry <= db->last_entry; i_entry++) {
    align_entry = db->align_array + i_entry;
    align_entry->orig_qual = db->orig_qual_area + i;
    align_entry->adj_qual = db->adj_qual_area + i;
    i += get_seq_length(i_entry) + incr;
  }
  if (i != t_size) fatalError("qual_area size");
}

reset_file(file)
     File *file;
{
  FILE *fp;
  int c;

  fp = file->fp;
  rewind(fp);

  if (!file->type) {
    while (isspace(c = fgetc(fp)));
    if (c != '>') {
      fprintf(stderr, "\nERROR: file %s not in FASTA format -- leading character not >",
	      file->name);
      exit(1);
    }
  }
  else { /* CALF file */
    while ((c = fgetc(file->fp)) && c != EOF); /* find end of ASCII section */
    if (c || fgetc(file->fp)) {
      fprintf(stderr, "\nERROR: file %s is not a .calf file of unaligned reads: ASCII section not followed by two 0 bytes", 
	      file->name);
      exit(1);
    }
  }
}

/* NB: following assumes FASTA format & that on entry fp is always pointing at the
   character directly behind the '>' at the beginning of the entry. On
   return, points to the appropriate
   buffer; thus to retain information when next entry is accessed, space
   must be allocated in the calling program. */

/* could speed up by avoiding rereading, except when necessary,; and by reading a line at a time */

get_next_file_entry(db)
     Database *db;
{
  int i, c, c_prev;
  unsigned char d;
  unsigned char *make_new_buffer();
  FILE *fp;
  char *id_buffer, *descrip_buffer;
  unsigned char *seq_buffer;
  int id_buffer_size, descrip_buffer_size, seq_buffer_size;

  fp = db->file->fp;
  id_buffer = db->id_buffer;
  descrip_buffer = db->descrip_buffer;
  seq_buffer = db->seq_buffer;
  id_buffer_size = db->id_buffer_size;
  descrip_buffer_size = db->descrip_buffer_size;
  seq_buffer_size = db->seq_buffer_size;

  id_pos = ftell(fp); 

  for (i = 0; !isspace(c = fgetc(fp)) && c != EOF; i++) {
    if (i < id_buffer_size) id_buffer[i] = c;
  }

  if (i >= id_buffer_size) {
    notify("id_buffer expanded\n");
    db->id_buffer = id_buffer = (char *)make_new_buffer((unsigned char *)id_buffer, i, 0);
    db->id_buffer_size = id_buffer_size = i;
    fseek(fp, id_pos, 0);
    for (i = 0; !isspace(c = fgetc(fp)) && c != EOF; i++) id_buffer[i] = c; 
  }

  id_buffer[i] = 0;

  if (c == EOF) return 0;

  for (; isspace(c) && c != '\n'; c = fgetc(fp));

  descrip_pos = ftell(fp) - 1;

  c_prev = c;
  for (i = 0; c != '\n' && c != EOF; c = fgetc(fp), i++) {
    if (i < descrip_buffer_size) descrip_buffer[i] = c; 
  }
  if (i >= descrip_buffer_size)  {
    notify("descrip_buffer expanded\n");
    db->descrip_buffer = descrip_buffer = (char *)make_new_buffer((unsigned char *)descrip_buffer, i, 0);
    db->descrip_buffer_size = descrip_buffer_size = i;
    fseek(fp, descrip_pos + 1, 0);
    c = c_prev;
    for (i = 0; c != '\n' && c != EOF; c = fgetc(fp), i++) descrip_buffer[i] = c; 
  }

  descrip_buffer[i] = 0;

  seq_pos = ftell(fp);

  for (i = 0; (c = fgetc(fp)) != '>' && c != EOF; ) {
    if (d = buffer_alph_mat[c]) {
      if (i < seq_buffer_size) seq_buffer[i] = d;

      i++;
    }
  }

  if (i >= seq_buffer_size)  {
    notify("seq_buffer expanded\n");
    db->seq_buffer = seq_buffer = make_new_buffer(seq_buffer, i + 1000, 1);
    db->seq_buffer_size = seq_buffer_size = i + 1000;
    fseek(fp, seq_pos, 0);
    for (i = 0; (c = fgetc(fp)) != '>' && c != EOF; ) if (d = buffer_alph_mat[c]) seq_buffer[i++] = d;
  }

  seq_buffer[i] = 0;

  db->seq_buffer_entry = db->comp_seq_buffer_entry = -1; /* SET THESE TO CORRECT VALUES !!! */ 

  db->seq_length = seq_length = i;

  return 1;
}

unsigned char *make_new_buffer(buffer, size, extend_flag)
     unsigned char *buffer;
     int size, extend_flag;
{
  int i;
  unsigned char *new_buffer;
  char *our_alloc();

  if (buffer)
    our_free(buffer - extend_flag);
  new_buffer = extend_flag + (unsigned char *)our_alloc((size + extend_flag + 1) * sizeof(unsigned char));
  new_buffer[-extend_flag] = 0;
  return new_buffer;
}

/*
char *make_new_buffer(buffer, add_size)
     char *buffer;
     int *add_size;
{
  int i;
  char *new_buffer;
  char *our_alloc();

  new_buffer = (char *)our_alloc((*add_size * 2 + 1) * sizeof(char));
  for (i = 0; i < *add_size; i++) new_buffer[i] = buffer[i];

  our_free(buffer);
  *add_size *= 2;
  return new_buffer;
}
*/

extern int col_convert[], n_c_cols;

/* initial read of database files to determine size needed, allocate space 
 store_comps is flag indicating whether complements are needed */
alloc_mem(db, store_comps)
     Database *db;
     int store_comps;
{
  FILE *fp;
  char *our_alloc();
  int c, n, n_entries, n_headers, dna_flag, nuc, q, flag;
  long int n_residues, comp_size, n_id_chars, n_descrip_chars;
  long int residues[256], uc_residues[256], col_convert_counts[256];

  n_id_chars = n_descrip_chars = 0;

  for (c = 0; c < 256; c++) residues[c] = uc_residues[c] = col_convert_counts[c] = 0;
  n_entries = 0;
  fp = db->file->fp;
  if (!db->file->type) {
    do {
      n_entries++;
      while (!isspace(c = fgetc(fp)) && c != EOF) n_id_chars++;
      for (; isspace(c) && c != '\n'; c = fgetc(fp));
      for (; c != '\n' && c != EOF; c = fgetc(fp)) n_descrip_chars++;
      while ((c = fgetc(fp)) != '>' && c != EOF) residues[c] += 1;
    } while (c != EOF);
  }
  else {
    n_headers = 0;
    for (c = fgetc(fp); c != EOF; c = fgetc(fp)) {
      if (c == 62) { /* start byte -- implies header */
	n_headers++;
	if (fgetc(fp)) fatalError("unaligned read header error: missing 0 byte");
	while (!isspace(c = fgetc(fp)) && c && c != EOF) n_id_chars++;
	for (; isspace(c); c = fgetc(fp));
	for (; c && c != EOF; c = fgetc(fp)) n_descrip_chars++;
	if (c) fatalError("premature EOF in unaligned reads");
	if (fgetc(fp) != 62) fatalError("unaligned read header error: missing start byte copy");
      }      
      else residues[c] += 1;
    }
    n_entries = residues[0]; /* # read ends */
  }

  printf("\n\nSequence file: %s    %d entries", db->file->name, n_entries);

  if (!db->file->type) {
    printf("\nResidue counts,  score matrix column # (* = last column):");

    for (c = n_residues = 0; c < 256; c++) {
      if (area_alph_mat[c] && residues[c]) {
	uc_residues[toupper(c)] += residues[c];
	n_residues += residues[c];
	col_convert_counts[col_convert[c]] += residues[c];
	printf("\n  %c    %8ld  %3d", (char)c, residues[c], col_convert[c] + 1);
	if (col_convert[c] == n_c_cols - 1)
	  printf("*");
      }
      else if (!isspace(c) && residues[c]) {
	printf("\n  %c    %8ld (presumed non-residue character -- ignored)", (char)c, residues[c]);
      }
    }
  }

  else {
    printf(" (%d with ASCII headers)\nResidue counts", n_headers);

    flag = 0;
    for (c = 1, n_residues = 0; c < 256; c++) {
      if (residues[c]) {
	q = c & 63;
	nuc = q && q < 62 ? "ACGT"[c >> 6] : c == 64 ? 'N' : c == 128 ? '-' : 'X';
	if (nuc == 'X') {
	  fprintf(stderr, "ERROR: %ld .calf read bytes with value %d", residues[c], c);
	  flag = 1;
	}
	uc_residues[nuc] += residues[c];
	n_residues += residues[c];
      }
    }
    for (c = 0; c < 256; c++) {
      if (uc_residues[c]) {
	printf("\n  %c    %8ld", (char)c, uc_residues[c]);
      }
    }
    if (flag)
      fatalError("anomalous .calf read bytes");
  }
  printf("\nTotal  %5ld residues", n_residues);

  for (c = n = 0; c < 256; c++) 
    if (isupper(c) && uc_residues[c] > .01 * n_residues) n++;

  printf("\n\n%d distinct alphabetic chars have freq > 1%% -- ", n);

  if (!db->file->type) {
    if (col_convert_counts[n_c_cols - 1] > .5 * n_residues
	|| !parameters->DNA_flag && n < 15
	|| parameters->DNA_flag == 1 && n > 6) {
      fprintf(stderr, "\n\nWARNING: POSSIBLE SCORE MATRIX/QUERY SEQUENCE INCOMPATIBILITY");
      fprintf(stdout, "\n\nWARNING: POSSIBLE SCORE MATRIX/QUERY SEQUENCE INCOMPATIBILITY");
    }
  }
    
  reset_file(db->file);
  /*
  rewind(fp);
  fgetc(fp);
  */

  db->t_length = n_residues;

  db->id_area_size = n_id_chars + n_entries + 1;
  size_conflict(db->id_area_size, n_id_chars + n_entries + 1, "id");
  db->id_area = (char *)our_alloc(db->id_area_size * sizeof(char));

  db->descrip_area_size = n_descrip_chars + n_entries + 1;
  size_conflict(db->descrip_area_size, n_descrip_chars + n_entries + 1, "descrip");
  db->descrip_area = (char *)our_alloc(db->descrip_area_size * sizeof(char));

  comp_size = n_residues + n_entries;
  if (store_comps) {
/* allow space for complements of entries in first file */
    n_entries *= 2;
    n_residues *= 2;
  }

  db->seq_area_size = n_residues + n_entries + 1;
  size_conflict(db->seq_area_size, n_residues + n_entries + 1, "seq");
  db->seq_area = (unsigned char *)our_alloc(db->seq_area_size * sizeof(unsigned char));

  db->comp_area = db->seq_area + comp_size;
  printf("\n\nAllocated space: %lu seqs, %lu ids, %lu descrips", 
	 (unsigned long int)db->seq_area_size, (unsigned long int)db->id_area_size, (unsigned long int)db->descrip_area_size);
}

size_conflict(size1, size2, label)
     SEQ_AREA size1;
     long int size2;
     char *label;
{
  SEQ_AREA sa_test;

  if (size1 == size2) return;
  printf("\nSize conflict with allocated %s data: %ld", label, size2);
  fprintf(stderr, "\nSize conflict with allocated %s data: %ld", label, size2);
  sa_test = -1;
  printf("\nVariable sizes: int: %d, long int: %d, SEQ_AREA: %d, signed: %d\n", (int)sizeof(int), (int)sizeof(long int), (int)sizeof(SEQ_AREA), sa_test < 0);
  fprintf(stderr, "\nVariable sizes: int: %d, long int: %d, SEQ_AREA: %d, signed: %d\n", 
	  (int)sizeof(int), (int)sizeof(long int), (int)sizeof(SEQ_AREA), sa_test < 0);
  fatalError("in query data allocation; may be fixable by changing SEQ_AREA typedef in swat.h and recompiling (see phrap.doc).");
}

file_to_mem(db) 
     Database *db;
{
  char *our_alloc();
  FILE *fp;
  int c, len, max_length;
  SEQ_AREA i_id, i_descrip, i_seq;
  unsigned char d;
  char *temp_ptr;
  unsigned char *seq_area; /* array containing all residue sequences; first char is 0, 
		     which is followed by each sequence in turn, each being 
		     terminated by 0 */
  char *id_area, *descrip_area; /* arrays containing ids and descriptions */

  seq_area = db->seq_area;
  id_area = db->id_area;
  descrip_area = db->descrip_area;

  id_area[0] = descrip_area[0] = seq_area[0] = 0;
  i_id = i_descrip = i_seq = 1;

  max_length = 0;
  
  fp = db->file->fp;
  if (!db->file->type) {
    do { 
      id_pos = i_id;
      while (!isspace(c = fgetc(fp)) && c != EOF) id_area[i_id++] = c;
      id_area[i_id++] = 0;
    
      for (; isspace(c) && c != '\n'; c = fgetc(fp));
    
      descrip_pos = i_descrip;
      for (; c != '\n' && c != EOF; c = fgetc(fp)) descrip_area[i_descrip++] = c;
      descrip_area[i_descrip++] = 0;
    
      seq_pos = i_seq;
    
      while ((c = fgetc(fp)) != '>' && c != EOF) 
	if (d = area_alph_mat[c]) seq_area[i_seq++] = d; 
    
      seq_length = i_seq - seq_pos;
      seq_area[i_seq++] = 0;
    
      append_seq_entry(db);
    
      if (seq_length > max_length) max_length = seq_length;
    
      temp_ptr = id_area + id_pos;
      len = strlen(temp_ptr);
      if (len >= 4 && !strcmp(temp_ptr + len - 4, ".seq")) {
	len -= 4;
	temp_ptr[len] = 0;
      }
    } while (c != EOF);
  }
  else { /* calf file; also get quality histogram (since no .qual file) */
    for (c = fgetc(fp); c != EOF; c = fgetc(fp)) {
      id_pos = i_id;
      descrip_pos = i_descrip;
      if (c == 62) { /* start byte -- implies header */
	fgetc(fp); /* 0 byte */
	while (!isspace(c = fgetc(fp)) && c) id_area[i_id++] = c;
	for (; isspace(c); c = fgetc(fp));
	for (; c; c = fgetc(fp)) descrip_area[i_descrip++] = c;
	fgetc(fp); /* start byte */
	c = fgetc(fp); /* 1st sequence byte */
      }      
      id_area[i_id++] = 0;
      descrip_area[i_descrip++] = 0;

      seq_pos = i_seq;
      for (; c; c = fgetc(fp)) 
	seq_area[i_seq++] = c; 

      seq_length = i_seq - seq_pos;
      seq_area[i_seq++] = 0;
    
      append_seq_entry(db);
    
      if (seq_length > max_length) max_length = seq_length;
    }
  }

  reset_file(db->file);
  
  if (i_id != db->id_area_size) {
    fprintf(stderr, "\n%lu %lu\n", (unsigned long int)i_id, (unsigned long int)db->id_area_size);
    fatalError("id_area_size");
  }
  if (i_descrip != db->descrip_area_size) fatalError("descrip_area_size");
  if (i_seq != (db->complements ? 1 + db->seq_area_size / 2 : db->seq_area_size)) {
    fprintf(stderr, "\n%d %lu %lu",db->complements, (unsigned long int)i_seq,  (unsigned long int)db->seq_area_size);
    fatalError("seq_area_size");
  }
/* make sure buffer is adequate to store complements */

/*
  if (max_length > db->seq_buffer_size) { 
    db->seq_buffer_size = max_length;
    our_free(db->seq_buffer);
    db->seq_buffer = (char *)our_alloc((db->seq_buffer_size + 1) * sizeof(char));
    db->seq_buffer[0] = 0;
    db->seq_buffer_entry = db->comp_seq_buffer_entry = -1;
  }
*/

/*  fprintf(stderr,"\n%s", seq_area + 1); */

/* following test probably not in correct location -- & NOW OBSOLETE */
  if (!strcmp(parameters->calling_program, "phrap") && max_length > 64000 && sizeof(Read_position) < 4)
    fatalError("Longest read sequence exceeds current limit (64000)\n -- use program phrap.longreads instead (create using make manyreads)");
}

set_score(entry_num, score) 
     int entry_num, score;
{
  Seq_entry *get_seq_entry();

  get_seq_entry(entry_num)->score = score;
}
 
Seq_entry *get_seq_entry(entry_num)
     int entry_num;
{
  int n;

  /*  if (entry_num < 0 || entry_num >= t_num_entries) fatalError("out of block range"); */
  n = entry_num / BLOCK_SIZE;
  return seq_block_table[n] + entry_num - n * BLOCK_SIZE; /* entry_num % BLOCK_SIZE; */
}

Database *get_db(entry_num)
     int entry_num;
{
  Database *db;

  for (db = head_db; db->first_entry > entry_num || db->first_entry < 0; db = db->next);
  return db; 
}

/* N.B. FOLLOWING ASSUME THAT id_buffer, descrip_buffer, AND seq_buffer ARE ALL LARGE ENOUGH */
char *get_id(entry_num)
     int entry_num;
{
  Seq_entry *seq_entry;
  Seq_entry *get_seq_entry();
  Database *db;
  Database *get_db();
  FILE *fp;
  int j, c;
  char *id_buffer;

  seq_entry = get_seq_entry(entry_num);
  db = get_db(entry_num);

  if (db->in_memory) return db->id_area + seq_entry->id_pos;

  fp = db->file->fp;
  fseek(fp, seq_entry->id_pos, 0);
  id_buffer = db->id_buffer;
  for (j = 0; !isspace(c = fgetc(fp)) && c != EOF;) id_buffer[j++] = c;
  id_buffer[j] = 0;
  return id_buffer;
}

char *get_descrip(entry_num)
     int entry_num;
{
  Seq_entry *seq_entry;
  Seq_entry *get_seq_entry();
  Database *db;
  Database *get_db();
  FILE *fp;
  int j, c;
  char *descrip_buffer;

  seq_entry = get_seq_entry(entry_num);
  db = get_db(entry_num);

  if (db->in_memory) return db->descrip_area + seq_entry->descrip_pos;

  fp = db->file->fp;
  fseek(fp, seq_entry->descrip_pos, 0);
  descrip_buffer = db->descrip_buffer;
  for (j = 0; (c = fgetc(fp)) != '\n' && c != EOF; ) descrip_buffer[j++] = c;
  descrip_buffer[j] = 0;
  return descrip_buffer;
}

unsigned char *get_seq(entry_num)
     int entry_num;
{
  Seq_entry *seq_entry;
  Seq_entry *get_seq_entry();
  Database *db;
  Database *get_db();
  FILE *fp;
  int j, c, length;
  POS seq_pos;
  unsigned char *seq_buffer;
  unsigned char *make_new_buffer();

  db = get_db(entry_num);
  seq_buffer = db->seq_buffer;
  if (db->seq_buffer_entry == entry_num) 
    return seq_buffer;

  seq_entry = get_seq_entry(entry_num);
  length = seq_entry->seq_length;
  seq_pos = seq_entry->seq_pos;

  if (db->in_memory) 
    return db->seq_area + seq_pos;

  if (db->seq_buffer_size < length) {
    notify("seq_buffer expanded\n");

    seq_buffer = db->seq_buffer = make_new_buffer(seq_buffer, length + 1000, 1);
    db->seq_buffer_size = length + 1000;
  }

  fp = db->file->fp;
  fseek(fp, seq_pos, 0);
  for (j = 0; j < length; ) {
    if (c = buffer_alph_mat[fgetc(fp)]) seq_buffer[j++] = c;
  }
  seq_buffer[j] = 0;
  db->seq_buffer_entry = entry_num;
  db->seq_length = length; /* IS THIS CORRECT */
  return seq_buffer;
}

unsigned char *get_comp_seq(entry_num)
     int entry_num;
{
  Seq_entry *seq_entry;
  Seq_entry *get_seq_entry();
  Database *db;
  Database *get_db();
  FILE *fp;
  int j, c, length;
  POS seq_pos;
  unsigned char *seq;
  unsigned char *comp_seq_buffer;
  unsigned char *make_new_buffer();

  db = get_db(entry_num);
  comp_seq_buffer = db->comp_seq_buffer;

  if (db->comp_seq_buffer_entry == entry_num) 
    return comp_seq_buffer;

  seq_entry = get_seq_entry(entry_num);
  length = seq_entry->seq_length;
  seq_pos = seq_entry->seq_pos;

  if (db->in_memory && db->complements)
    return db->seq_area + (db->seq_area_size - seq_pos - length);

  if (db->comp_seq_buffer_size < length) {
    notify("comp_seq_buffer expanded\n");
    comp_seq_buffer = db->comp_seq_buffer = make_new_buffer(comp_seq_buffer, length + 1000, 1);
    db->comp_seq_buffer_size = length + 1000;
  }

  comp_seq_buffer[length] = 0;

  if (db->in_memory) {
    seq = db->seq_area + seq_pos;
    for (j = length - 1; j >= 0; j--, seq++) 
      comp_seq_buffer[j] = buffer_comp_mat[*seq];
  }
  else { /* read from file into comp_seq_buffer */
    fp = db->file->fp;
    fseek(fp, seq_pos, 0);
    for (j = length - 1; j >= 0; ) {
      if (c = buffer_alph_mat[fgetc(fp)]) comp_seq_buffer[j--] = buffer_comp_mat[c];
    }
    /* notify("COMP_SEQ_BUFFER FILLED"); */
  }

  db->comp_seq_length = length;
  db->comp_seq_buffer_entry = entry_num;
  return comp_seq_buffer;
}

int get_seq_length(entry_num)
     int entry_num;
{
  Seq_entry *get_seq_entry();

  return get_seq_entry(entry_num)->seq_length;
}

Cand_pair *get_cand_pairs(entry_num)
     int entry_num;
{
  Seq_entry *get_seq_entry();

  return get_seq_entry(entry_num)->cand_pairs;
}

Aligned_pair *get_aligned_pairs(entry_num)
     int entry_num;
{
  Seq_entry *get_seq_entry();

  return get_seq_entry(entry_num)->aligned_pairs;
}

Align_info *get_align_entry(entry_num)
     int entry_num;
{
  Database *get_db();

  return get_db(entry_num)->align_array + entry_num;
}

char *get_orig_qual(entry_num)
     int entry_num;
{
  Database *get_db();

  return get_db(entry_num)->align_array[entry_num].orig_qual;
}

int get_equiv_class(entry_num)
     int entry_num;
{
  Database *get_db();

  return get_db(entry_num)->align_array[entry_num].equiv_class;
}

set_equiv_class(entry_num, class)
     int entry_num, class;
{
  Database *get_db();

  get_db(entry_num)->align_array[entry_num].equiv_class = class;
}

char *get_adj_qual(entry_num)
     int entry_num;
{
  Database *get_db();

  return get_db(entry_num)->align_array[entry_num].adj_qual;
}

get_chemistry(entry_num)
     int entry_num;
{
  Database *get_db();

  return (int)(get_db(entry_num)->align_array[entry_num].chemistry);
}

get_read_direction(entry_num)
     int entry_num;
{
  Database *get_db();

  return (int)(get_db(entry_num)->align_array[entry_num].direction);
}

set_chemistry(entry_num, value)
     int entry_num, value;
{
  Database *get_db();

  get_db(entry_num)->align_array[entry_num].chemistry = value;
}

set_direction(entry_num, value)
     int entry_num, value;
{
  Database *get_db();

  get_db(entry_num)->align_array[entry_num].direction = value;
}

int append_seq_entry(db)
     Database *db;
{
  Seq_entry *seq_entry;
  int i_entry;

  i_entry = t_num_entries % BLOCK_SIZE;

  if (!i_entry) alloc_seq_block();

  seq_entry = seq_block_table[num_blocks - 1] + i_entry;
  seq_entry->id_pos = id_pos;
  seq_entry->descrip_pos = descrip_pos;
  seq_entry->seq_pos = seq_pos;
  seq_entry->seq_length = seq_length;
  seq_entry->aligned_pairs = 0;
  seq_entry->cand_pairs = 0;
  seq_entry->query_domains = 0;
  seq_entry->score = 0;
  if (db->first_entry < 0) db->first_entry = t_num_entries;
  db->last_entry = t_num_entries;
  db->num_entries += 1;
  t_num_entries++;
  return t_num_entries - 1;
}

remove_seq_entry(db)
     Database *db;
{
  int i_entry;

  db->last_entry -= 1;
  db->num_entries -= 1;

  t_num_entries--;

  i_entry = t_num_entries % BLOCK_SIZE;

  /* NEED TO REMOVE TAGS ALSO! */
  if (!i_entry) remove_seq_block();
}

alloc_seq_block()
{
  char *our_alloc();

  if (num_blocks >= MAX_NUM_BLOCKS) fatalError("MAX_NUM_BLOCKS exceeded");
  if (sizeof(Entry_index) < 4 && (num_blocks + 1) * BLOCK_SIZE > 64000) /* NOW OBSOLETE */
    fatalError("Number of sequences exceeds current limit (64000)\n -- use .manyreads version instead (create using make manyreads)");

  seq_block_table[num_blocks++] = (Seq_entry *)our_alloc(BLOCK_SIZE * sizeof(Seq_entry));
  /* notify("BLOCK ALLOCATED"); */
}

remove_seq_block()
{
  num_blocks--;
  our_free(seq_block_table[num_blocks]);
  seq_block_table[num_blocks] = 0;
}

write_entry(fp, id, descrip, seq)
     FILE *fp;
     char *id, *descrip;
     unsigned char *seq;
{
  int j;

  fprintf(fp,">%s", id);
  if (descrip) fprintf(fp,"  %s", descrip);

  if (parameters->DNA_flag >= 3) {
    for (j = 0; ; j++) {
      if (!(j % 50)) fprintf(fp,"\n");
      if (!seq[j]) break;
      fprintf(fp,"%c", seq[j] & 63 ? "ACGT"[seq[j] >> 6] : " N-*"[seq[j] >> 6]);
    }
  }
  else {
    for (j = 0; ; j++) {
      if (!(j % 50)) fprintf(fp,"\n");
      if (!seq[j]) break;
      fprintf(fp,"%c", seq[j]);
    }
  }

  fprintf(fp,"\n");
}
