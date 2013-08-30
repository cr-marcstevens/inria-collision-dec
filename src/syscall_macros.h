/**
 * \file syscall_macros.h
 * \brief Syscall wrapper
 *
 * Same syntax as original syscall but uppercased. Abort with a nice error message on failure.
 *
 * \author Grégory Landais <gregory.landais@inria.fr>
 */
#ifndef SYSCALL_MACROS_H
#define SYSCALL_MACROS_H

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>

#define xerror(file, line, msg) do{fprintf(stderr, "%s:%d %s : %s\n", file, line, (msg), strerror(errno)); exit(EXIT_FAILURE);}while(0)

static inline void* xmalloc(size_t size, int line, char* file) {
	void* ptr;
	ptr = malloc(size);
	if(!ptr)
		xerror(file, line, "malloc");
	return ptr;
}

static inline void* xcalloc(size_t nmemb, size_t size, int line, char* file) {
	void* ptr;
	ptr = calloc(nmemb, size);
	if(!ptr)
		xerror(file, line, "calloc");
	return ptr;
}

static inline void* xrealloc(void *ptr, size_t size, int line, char* file) {
	void* new;
	new = realloc(ptr, size);
	if (!new && size != 0)
		xerror(file, line, "realloc");
	return new;
}

static inline FILE *xfopen(const char *filename, char *mode, int line, char* file) {
	FILE *f;
	f = fopen(filename, mode);
	if (!f)
		xerror(file, line, "fopen");
	return f;
}

/* not useful and not really wanted */
/*
static inline int xfread(void *ptr, size_t size, size_t nmemb, FILE *f, int line, char* file) {
	size_t n;
	n = fread(ptr, size, nmemb, f);
	if(n != nmemb)
		xerror(file, line, "fread");
	return n;
}

static inline void xfwrite(void *ptr, size_t size, size_t nmemb, FILE *f, int line, char* file) {
	size_t n;
	n = fwrite(ptr, size, nmemb, f);
	if (n != nmemb)
		xerror(file, line, "fwrite");
}

static inline void xfclose(FILE *f, int line, char* file) {
	if (fclose(f) == EOF)
		xerror(file, line, "fclose");
}
*/
#define MALLOC(n) xmalloc((n), __LINE__, __FILE__)
#define CALLOC(n, m) xcalloc((n), (m), __LINE__, __FILE__)
#define REALLOC(p, n) xrealloc((p), (n), __LINE__, __FILE__)
#define FOPEN(f, m) xfopen((f), (m), __LINE__, __FILE__)
/* for consistency */
#define FREAD(p, n, m, f) fread((p), (n), (m), (f))
#define FWRITE(p, n, m, f) fwrite((p), (n), (m), (f))
#define FCLOSE(f) fclose((f))

#endif
