extern int GetZ(dictionary *scan_file_dict);
extern int GetZplusOrZminus(dictionary *scan_file_dict, const char *sectionname, const char *nf_or_ff, char *result);
extern int CreateZlisting(dictionary *scan_file_dict, const char *sectionname_z1, const char *sectionname_z2, const char *ZplusOrZminus);
