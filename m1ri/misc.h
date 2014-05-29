#ifndef MISC_H
#define MISC_H

/** *
 * \brief Print error message and abort(). 
 * 
 * The function accepts additional
 * parameters like printf, so e.g. m1ri_die("foo %d bar %f\n",1 ,2.0)
 * is valid and will print the string "foo 1 bar 2.0" before dying.
 *
 * \param errormessage a string to be printed.
 *
 * \warning The provided string is not free'd.
 */

void m1ri_die(const char *errormessage, ...);


#endif //MISC_H
