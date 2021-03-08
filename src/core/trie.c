/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2003-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge, MA 02139 USA

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include "igraph_types.h"
#include "igraph_memory.h"
#include "igraph_random.h"
#include "igraph_error.h"

#include "core/trie.h"

#include "config.h"

#include <string.h>

/**
 * \ingroup igraphtrie
 * \brief Creates a trie node (not to be called directly)
 * \return Error code: errors by igraph_strvector_init(),
 *         igraph_vector_ptr_init() and igraph_vector_init() might be returned.
 */

static int igraph_i_trie_init_node(igraph_trie_node_t *t) {
    IGRAPH_STRVECTOR_INIT_FINALLY(&t->strs, 0);
    IGRAPH_VECTOR_PTR_INIT_FINALLY(&t->children, 0);
    IGRAPH_VECTOR_INIT_FINALLY(&t->values, 0);
    IGRAPH_FINALLY_CLEAN(3);
    return 0;
}

static void igraph_i_trie_destroy_node(igraph_trie_node_t *t);

/**
 * \ingroup igraphtrie
 * \brief Creates a trie.
 * \return Error code: errors by igraph_strvector_init(),
 *         igraph_vector_ptr_init() and igraph_vector_init() might be returned.
 */

int igraph_trie_init(igraph_trie_t *t, igraph_bool_t storekeys) {
    t->maxvalue = -1;
    t->storekeys = storekeys;
    IGRAPH_CHECK(igraph_i_trie_init_node( (igraph_trie_node_t *) t ));
    IGRAPH_FINALLY(igraph_i_trie_destroy_node, (igraph_trie_node_t *) t );
    if (storekeys) {
        IGRAPH_CHECK(igraph_strvector_init(&t->keys, 0));
    }

    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}

/**
 * \ingroup igraphtrie
 * \brief Destroys a node of a trie (not to be called directly).
 */

static void igraph_i_trie_destroy_node_helper(igraph_trie_node_t *t, igraph_bool_t sfree) {
    long int i;
    igraph_strvector_destroy(&t->strs);
    for (i = 0; i < igraph_vector_ptr_size(&t->children); i++) {
        igraph_trie_node_t *child = VECTOR(t->children)[i];
        if (child != 0) {
            igraph_i_trie_destroy_node_helper(child, 1);
        }
    }
    igraph_vector_ptr_destroy(&t->children);
    igraph_vector_destroy(&t->values);
    if (sfree) {
        IGRAPH_FREE(t);
    }
}

static void igraph_i_trie_destroy_node(igraph_trie_node_t *t) {
    igraph_i_trie_destroy_node_helper(t, 0);
}

/**
 * \ingroup igraphtrie
 * \brief Destroys a trie (frees allocated memory).
 */

void igraph_trie_destroy(igraph_trie_t *t) {
    if (t->storekeys) {
        igraph_strvector_destroy(&t->keys);
    }
    igraph_i_trie_destroy_node( (igraph_trie_node_t*) t);
}


/**
 * \ingroup igraphtrie
 * \brief Internal helping function for igraph_trie_t
 */

static long int igraph_i_strdiff(const char *str, const char *key) {

    long int diff = 0;
    while (key[diff] != '\0' && str[diff] != '\0' && str[diff] == key[diff]) {
        diff++;
    }
    return diff;
}

/**
 * \ingroup igraphtrie
 * \brief Search/insert in a trie (not to be called directly).
 *
 * @return Error code:
 *         - <b>IGRAPH_ENOMEM</b>: out of memory
 */

int igraph_trie_get_node(igraph_trie_node_t *t, const char *key,
                         igraph_real_t newvalue, long int *id) {
    char *str;
    long int i;
    igraph_bool_t add;

    /* If newvalue is negative, we don't add the node if nonexistent, only check
     * for its existence */
    add = (newvalue >= 0);

    for (i = 0; i < igraph_strvector_size(&t->strs); i++) {
        long int diff;
        igraph_strvector_get(&t->strs, i, &str);
        diff = igraph_i_strdiff(str, key);

        if (diff == 0) {

            /* ------------------------------------ */
            /* No match, next */

        } else if (str[diff] == '\0' && key[diff] == '\0') {

            /* ------------------------------------ */
            /* They are exactly the same */
            if (VECTOR(t->values)[i] != -1) {
                *id = (long int) VECTOR(t->values)[i];
                return 0;
            } else {
                VECTOR(t->values)[i] = newvalue;
                *id = (long int) newvalue;
                return 0;
            }

        } else if (str[diff] == '\0') {

            /* ------------------------------------ */
            /* str is prefix of key, follow its link if there is one */
            igraph_trie_node_t *node = VECTOR(t->children)[i];
            if (node != 0) {
                return igraph_trie_get_node(node, key + diff, newvalue, id);
            } else if (add) {
                igraph_trie_node_t *node = IGRAPH_CALLOC(1, igraph_trie_node_t);
                if (node == 0) {
                    IGRAPH_ERROR("cannot add to trie", IGRAPH_ENOMEM);
                }
                IGRAPH_STRVECTOR_INIT_FINALLY(&node->strs, 1);
                IGRAPH_VECTOR_PTR_INIT_FINALLY(&node->children, 1);
                IGRAPH_VECTOR_INIT_FINALLY(&node->values, 1);
                IGRAPH_CHECK(igraph_strvector_set(&node->strs, 0, key + diff));
                VECTOR(node->children)[0] = 0;
                VECTOR(node->values)[0] = newvalue;

                VECTOR(t->children)[i] = node;

                *id = (long int) newvalue;
                IGRAPH_FINALLY_CLEAN(3);
                return 0;
            } else {
                *id = -1;
                return 0;
            }

        } else if (key[diff] == '\0' && add) {

            /* ------------------------------------ */
            /* key is prefix of str, the node has to be cut */
            char *str2;

            igraph_trie_node_t *node = IGRAPH_CALLOC(1, igraph_trie_node_t);
            if (node == 0) {
                IGRAPH_ERROR("cannot add to trie", IGRAPH_ENOMEM);
            }
            IGRAPH_STRVECTOR_INIT_FINALLY(&node->strs, 1);
            IGRAPH_VECTOR_PTR_INIT_FINALLY(&node->children, 1);
            IGRAPH_VECTOR_INIT_FINALLY(&node->values, 1);
            IGRAPH_CHECK(igraph_strvector_set(&node->strs, 0, str + diff));

            VECTOR(node->children)[0] = VECTOR(t->children)[i];
            VECTOR(node->values)[0] = VECTOR(t->values)[i];

            str2 = strdup(str);
            if (str2 == 0) {
                IGRAPH_ERROR("cannot add to trie", IGRAPH_ENOMEM);
            }
            str2[diff] = '\0';
            IGRAPH_FINALLY(igraph_free, str2);
            IGRAPH_CHECK(igraph_strvector_set(&t->strs, i, str2));
            IGRAPH_FREE(str2);
            IGRAPH_FINALLY_CLEAN(4);

            VECTOR(t->values)[i] = newvalue;
            VECTOR(t->children)[i] = node;

            *id = (long int) newvalue;
            return 0;

        } else if (add) {

            /* ------------------------------------ */
            /* the first diff characters match */
            char *str2;

            igraph_trie_node_t *node = IGRAPH_CALLOC(1, igraph_trie_node_t);
            if (node == 0) {
                IGRAPH_ERROR("cannot add to trie", IGRAPH_ENOMEM);
            }
            IGRAPH_STRVECTOR_INIT_FINALLY(&node->strs, 2);
            IGRAPH_VECTOR_PTR_INIT_FINALLY(&node->children, 2);
            IGRAPH_VECTOR_INIT_FINALLY(&node->values, 2);
            IGRAPH_CHECK(igraph_strvector_set(&node->strs, 0, str + diff));
            IGRAPH_CHECK(igraph_strvector_set(&node->strs, 1, key + diff));
            VECTOR(node->children)[0] = VECTOR(t->children)[i];
            VECTOR(node->children)[1] = 0;
            VECTOR(node->values)[0] = VECTOR(t->values)[i];
            VECTOR(node->values)[1] = newvalue;

            str2 = strdup(str);
            if (str2 == 0) {
                IGRAPH_ERROR("cannot add to trie", IGRAPH_ENOMEM);
            }
            str2[diff] = '\0';
            IGRAPH_FINALLY(igraph_free, str2);
            IGRAPH_CHECK(igraph_strvector_set(&t->strs, i, str2));
            IGRAPH_FREE(str2);
            IGRAPH_FINALLY_CLEAN(4);

            VECTOR(t->values)[i] = -1;
            VECTOR(t->children)[i] = node;

            *id = (long int) newvalue;
            return 0;
        } else {

            /* ------------------------------------------------- */
            /* No match, but we requested not to add the new key */
            *id = -1;
            return 0;
        }
    }

    /* ------------------------------------ */
    /* Nothing matches */

    if (add) {
        IGRAPH_CHECK(igraph_vector_ptr_reserve(&t->children,
                                               igraph_vector_ptr_size(&t->children) + 1));
        IGRAPH_CHECK(igraph_vector_reserve(&t->values, igraph_vector_size(&t->values) + 1));
        IGRAPH_CHECK(igraph_strvector_add(&t->strs, key));

        igraph_vector_ptr_push_back(&t->children, 0); /* allocated */
        igraph_vector_push_back(&t->values, newvalue); /* allocated */
        *id = (long int) newvalue;
    } else {
        *id = -1;
    }

    return 0;
}

/**
 * \ingroup igraphtrie
 * \brief Search/insert in a trie.
 */

int igraph_trie_get(igraph_trie_t *t, const char *key, long int *id) {
    if (!t->storekeys) {
        IGRAPH_CHECK(igraph_trie_get_node( (igraph_trie_node_t*) t,
                                           key, t->maxvalue + 1, id));
        if (*id > t->maxvalue) {
            t->maxvalue = *id;
        }
        return 0;
    } else {
        int ret;
        igraph_error_handler_t *oldhandler;
        oldhandler = igraph_set_error_handler(igraph_error_handler_ignore);
        /* Add it to the string vector first, we can undo this later */
        ret = igraph_strvector_add(&t->keys, key);
        if (ret != 0) {
            igraph_set_error_handler(oldhandler);
            IGRAPH_ERROR("cannot get element from trie", ret);
        }
        ret = igraph_trie_get_node( (igraph_trie_node_t*) t,
                                    key, t->maxvalue + 1, id);
        if (ret != 0) {
            igraph_strvector_resize(&t->keys, igraph_strvector_size(&t->keys) - 1);
            igraph_set_error_handler(oldhandler);
            IGRAPH_ERROR("cannot get element from trie", ret);
        }

        /* everything is fine */
        if (*id > t->maxvalue) {
            t->maxvalue = *id;
        } else {
            igraph_strvector_resize(&t->keys, igraph_strvector_size(&t->keys) - 1);
        }
        igraph_set_error_handler(oldhandler);
    }

    return 0;
}

/**
 * \ingroup igraphtrie
 * \brief Search/insert in a trie (for internal use).
 *
 * @return Error code:
 *         - <b>IGRAPH_ENOMEM</b>: out of memory
 */

int igraph_trie_get2(igraph_trie_t *t, const char *key, long int length,
                     long int *id) {
    char *tmp = IGRAPH_CALLOC(length + 1, char);

    if (tmp == 0) {
        IGRAPH_ERROR("Cannot get from trie", IGRAPH_ENOMEM);
    }

    strncpy(tmp, key, length);
    tmp[length] = '\0';
    IGRAPH_FINALLY(igraph_free, tmp);
    IGRAPH_CHECK(igraph_trie_get(t, tmp, id));
    IGRAPH_FREE(tmp);
    IGRAPH_FINALLY_CLEAN(1);
    return 0;
}

/**
 * \ingroup igraphtrie
 * \brief Search in a trie.
 * This variant does not add \c key to the trie if it does not exist.
 * In this case, a negative id is returned.
 */

int igraph_trie_check(igraph_trie_t *t, const char *key, long int *id) {
    IGRAPH_CHECK(igraph_trie_get_node( (igraph_trie_node_t*) t,
                                       key, -1, id));
    return 0;
}

/**
 * \ingroup igraphtrie
 * \brief Get an element of a trie based on its index.
 */

void igraph_trie_idx(igraph_trie_t *t, long int idx, char **str) {
    igraph_strvector_get(&t->keys, idx, str);
}

/**
 * \ingroup igraphtrie
 * \brief Returns the size of a trie.
 */

long int igraph_trie_size(igraph_trie_t *t) {
    return t->maxvalue + 1;
}

/* Hmmm, very dirty.... */

int igraph_trie_getkeys(igraph_trie_t *t, const igraph_strvector_t **strv) {
    *strv = &t->keys;
    return 0;
}
