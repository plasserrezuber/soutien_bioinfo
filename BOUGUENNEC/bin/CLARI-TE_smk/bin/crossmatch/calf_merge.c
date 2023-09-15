/*****************************************************************************
#   Copyright (C) 2008 by Phil Green.                          
#   All rights reserved.                           
#                                                                           
#   This software is part of a beta-test version of the swat/cross_match/phrap 
#   package.  It should not be redistributed, or
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

/*
    program to merge and/or modify calf files, and perform integrity checking. Merged
   file has (i) initial ASCII section which consists of the command
   line used to invoke this program, followed by the concatenated
   ASCII sections of the input files (separated by intervening \n's);
   (ii) the merged alignments; (iii) the concatenated unaligned read sections.

   All input files must have the same reference sequences, in the same
   order; the only exception being that some input files can be ASCII
   only, with no alignment or unaligned read sections (such files
   should have no terminating 0 byte).

   command line: calf_merge -out output_filename input_filename1 [input_filename2] [input_filename3] [input_filename4] ... [-option1] [-option2]

   (in fact the arguments to calf_merge can appear in any order, the only requirement being that the output_filename must immediately follow -out)

   Note that it is permissible to have a single input_file -- this is useful for performing integrity checks, stripping ASCII headers, etc.

    currently available options:
       -delete_type3 : causes all uncovered reference segments to be output as type 2 records 
          (size only), rather than type 3 (packed sequence). Note that one can do the reverse  
           (convert all size-only to packed refseq) by merging with a
           CALF file that has only the refseqs (as a series of type 2 records, one for each
           reference sequence) with no reads
       -delete_file_ASCII : deletes ASCII file sections. This is useful when changes
          to this section need to be made -- a revised version of the ASCII section can be merged 
          back in later by executing calf_merge with one of the files consisting entirely
          of an ASCII header (such a file should not contain a terminating 0 byte).
       -only_file_ASCII : deletes the alignments and unaligned reads (i.e. everything EXCEPT the ASCII file sections).
       -delete_read_ASCII : deletes ASCII read headers. 
       -delete_unaligned : delete all unaligned reads
       -delete_alignments : delete refseq and aligned reads 
       -delete_ref_offsets : do not include reference co-ord offsets in any continuation pointer in the output file
       -all_ref_offsets : include reference co-ord offsets in all continuation pointers in the output file

   if neither -delete_ref_offsets nor -all_ref_offsets is specified, then reference offsets
       in the input files are preserved

   all files are kept open during the run, so the total # must not exceed FOPEN_MAX. A check
   for this is performed.
 
   TO BE ADDED: 

   more error checks -- eg. 
       (i) check read operations for unexpected EOF
       (ii) check read operations for consistency of refseq locationss, and when some at new ref & others not (this not nec an error)    
        (iii) for fatal errors, output file_stat info

   add options to 
     (i) remove redundant reads
     (ii) when there is only one input file, shrink pointer bytes to min possible
     (iii) extract aligned reads (deleting reference sequence)
     (iv) merge a file containing unaligned reads only (can always just concatenate)

   need utilities to 
      (i) slide alignments in gapped regions after merge -- to deal with
        placements of nucs in long reference gaps (partic at beginnings of refseqs)

   co-ord system for ref seqs needs to be described in format description
*/

#include "calf.h"

extern File_stat *head_file_stat;

main(argc,argv)
     int argc;
     char *argv[];
{
  char output_file_name[101];
  FILE *fp_out;
  File_stat *file_stat;
  int i, gap_col, last_col, ref_offsets, n_files;
  int c, type, prev_type, ref_nuc, prev_nuc, half, uncov, new_ref;
  int delete_type3, delete_file_ASCII, only_file_ASCII, delete_read_ASCII, delete_unaligned, delete_alignments, keep;
  long int ref_loc, prev_loc;

  fp_out = 0;
  delete_type3 = delete_file_ASCII = only_file_ASCII = delete_read_ASCII = delete_unaligned = delete_alignments = 0;
  ref_offsets = 1;
  for (i = 1; i < argc; i++) {
    if (!strcmp("-out", argv[i])) {
      if (i >= argc - 1) fatal_error("no output file indicated");
      if (fp_out) fatal_error("more than one output file indicated");
      i++;
      fp_out = fopen(argv[i], "w+b"); /* open for read + write binary access */
      if (!fp_out) fatal_error("cannot open output file");
      strncpy(output_file_name, argv[i], 100);
    }
    else if (!strcmp("-delete_type3", argv[i])) {
      delete_type3 = 1;
    }
    else if (!strcmp("-delete_file_ASCII", argv[i])) {
      delete_file_ASCII = 1;
    }
    else if (!strcmp("-only_file_ASCII", argv[i])) {
      only_file_ASCII = 1;
    }
    else if (!strcmp("-delete_read_ASCII", argv[i])) {
      delete_read_ASCII = 1;
    }
    else if (!strcmp("-delete_unaligned", argv[i])) {
      delete_unaligned = 1;
    }
    else if (!strcmp("-delete_alignments", argv[i])) {
      delete_alignments = 1;
    }
    else if (!strcmp("-delete_ref_offsets", argv[i])) {
      ref_offsets = 0;
    }
    else if (!strcmp("-all_ref_offsets", argv[i])) {
      ref_offsets = 2;
    }
    else if (argv[i][0] == '-') {
      fatal_error("unrecognized option");
    }
    else {
      append_file_stat(argv[i]);
    }
  }
  if (!fp_out) fatal_error("no output file indicated"); /* maybe OK?? */

  for (file_stat = head_file_stat; file_stat; file_stat = file_stat->next)
    if (!strcmp(file_stat->name, output_file_name)) fatal_error("identically named input & output files");
 

  /* write initial ASCII line indicating command line used */

  for (i = 0; i < argc; i++) fprintf(fp_out, "%s ", argv[i]);

  /* read ASCII parts of files, & write to output (concatenated, with intervening \n's) */

  keep = !delete_file_ASCII;
  for (file_stat = head_file_stat, n_files = 0; file_stat; file_stat = file_stat->next) {
    if (keep) fprintf(fp_out, "\n");
    while ((c = fgetc(file_stat->fp)) && c != EOF) 
      if (keep) fputc(c, fp_out);
    if (c == EOF) {
      file_stat->ascii_only = 1;
    }
    else n_files++;
  }
  if (only_file_ASCII || !n_files) exit(1); /* finish with no terminating 0 */

  fputc(0, fp_out); /* terminates ASCII in output */

  /* now read the alignment portion of files, modify where necessary & write to output.

   strategy: for each input file, have current record, and refseq location (= co-ord of ref
     seq base (origin 1) for current record (the 1st base in the unaligned region, 
      in the case of a type 2 record); gaps in ref do not increment 
     for type 3 record, half (0 or 1) indicates whether in 1st half or 2d half of byte
     # active reads (kept by counting start/end bytes) -- nec when inserting gap

     in each round, output a record (or increment a type 3 record) unless in size-only uncovered region.
     if any input files have type 1 record at that pos'n
           terminate any active type 2 or type 3 record, 
               initiate type 1, 
               process type 1 records
               insert gap bytes for files with active reads but no records (arises when there
                 is a new gap)
                            
      if not, any type 3 record data at that pos'n -- 
                   append to or initiate growing type 3 record
  */

  keep = !delete_alignments;

  for (file_stat = head_file_stat; file_stat; file_stat = file_stat->next) {
    if (file_stat->ascii_only) continue;
    get_record_header(file_stat);
  }

  prev_type = prev_loc = 0;

  for (;;) {

    /* find minimum ref seq location of any current input record */
    ref_loc = (long)1 << 50; /* bigger than any that should occur */
    for (file_stat = head_file_stat, new_ref = 0; file_stat; file_stat = file_stat->next) {
      if (file_stat->ascii_only) continue;
      new_ref += file_stat->prev_type;
      if (ref_loc > file_stat->ref_loc /* always satisfied for first file_stat */
	  || ref_loc == file_stat->ref_loc
	  && type != file_stat->type
	  && (file_stat->type == 1 
	      || file_stat->type == 3 && type != 1
	      || file_stat->type == 2 && !type)) { /* preferred types: 1, over 3, over 2, over 0 */
	type = file_stat->type;
	ref_loc = file_stat->ref_loc;
	ref_nuc = file_stat->ref_nuc;
      } 
    }    

    new_ref = !new_ref; /* = 1 if all inputs at new ref seq */

    /* first finish up any uncovered region records */
    if (prev_type == 3 && (new_ref || type != 3)) { /* finish up type 3 record */
      if (keep) {
	if (half) fputc(prev_nuc << 4, fp_out); 
	fputc(0, fp_out); 
      }
    }

    gap_col = type == 1 && !ref_nuc; /* gap in reference */

    if (new_ref || type != 2) {
      last_col = gap_col ? ref_loc : ref_loc - 1; /* prev ref co-ord -- accounting for gap case */
      uncov = last_col - prev_loc;
      if (new_ref) uncov--; /* to cancel added base at beginning of new ref */

      if (uncov) { /* missing segment of reference; insert type 2 record */
	if (keep) {
	  fputc((prev_type << 2) + 2, fp_out); /* header byte */
	  write_offset(fp_out, 2, 0, uncov);
	  fputc(0, fp_out);
	}
	prev_type = 2;
	prev_loc = last_col; 
      }
    }

    if (new_ref) prev_type = 0;

    /* now deal with current records */

    if (type == 1) { /* write new type 1 record */
      fputc((ref_nuc << 4) + (prev_type << 2) + 1, fp_out); /* header byte */
      prev_type = 1;
      prev_loc = ref_loc;    
    }
    else if (type == 3) { /* initiate or extend current type 3 record */
      if (!delete_type3) { 
	if (prev_type != 3) {
	  if (keep) 
	    fputc((prev_type << 2) + 3, fp_out); /* header byte */
	  half = 0;
	  prev_type = 3;
	}
	if (half) {
	  if (keep) 
	    fputc((prev_nuc << 4) + ref_nuc, fp_out); 
	}
	prev_nuc = ref_nuc;
	prev_loc = ref_loc;    
	half = !half;
      }
    }
    else if (type == 2) { /* nothing need be done for type 2 */
      ;
    }
    else if (!type) { /* end of alignments; append empty record & break */
      for (file_stat = head_file_stat; file_stat; file_stat = file_stat->next) 
	if (!file_stat->ascii_only && file_stat->type)
	  fatal_error("premature alignment end");
      fputc(0, fp_out);
      break;
    }
    else fatal_error("unrecognized type");
    
    for (file_stat = head_file_stat; file_stat; file_stat = file_stat->next) {
      if (file_stat->ascii_only) continue;
      if (file_stat->ref_loc == ref_loc) {
	if (file_stat->type != 2 && file_stat->ref_nuc != ref_nuc) 
	  fatal_error("reference nucleotide incompatibility");
	if (file_stat->type == 1) read_write_type_1(file_stat, keep ? fp_out : (FILE *)0, ref_offsets, fp_out && !delete_read_ASCII);
	else if (file_stat->type == 2) read_type_2(file_stat);
	else if (file_stat->type == 3) read_type_3(file_stat, 0);
	else fatal_error("premature alignment end");
      }
      else if (file_stat->n_active) { /* unshared gap: write gap bytes */
	if (!gap_col) fatal_error("alignment parsing"); /* should only occur at gap */
	if (keep) 
	  for (i = 0; i < file_stat->n_active; i++)
	    fputc(128, fp_out); /* gap byte */
      }
    }
    if (type == 1) {
      if (keep) fputc(0, fp_out); /* terminate type 1 record */
    }
  } 

  if (delete_unaligned) exit(1);

  /* now write unaligned reads to output (concatenated) */
  for (file_stat = head_file_stat; file_stat; file_stat = file_stat->next) {
    if (file_stat->ascii_only) continue;

    if (delete_read_ASCII) {
      while ((c = fgetc(file_stat->fp)) != EOF) {
	if (c == 62) { /* start marker byte; must be followed by ASCII header */
	  if (fgetc(file_stat->fp)) fatal_error("missing 0 byte in unaligned read header");
	  do { 
	    c = fgetc(file_stat->fp);
	  } while (c); /* find next 0 byte, which ends header */
	  if (fgetc(file_stat->fp) != 62)
	    fatal_error("missing start byte in unaligned read header");
	} 
	else {
	  fputc(c, fp_out);
	}
      }
    }
    else {
      while ((c = fgetc(file_stat->fp)) != EOF) 
	fputc(c, fp_out);
    }

    fclose(file_stat->fp);
  }

  fclose(fp_out);
}


