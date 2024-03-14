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

#include "core/trie.h"
#include "internal/hacks.h" /* strdup */

#include <assert.h>
#include <string.h>


/*
 * igraph_trie_t is a data structures that stores an ordered list of strings.
 * It allows an efficient lookup of the index of a string. It has the capability
 * to also store the list of strings directly for reverse lookup of strings
 * by index.
 */

/* Allocates memory for a trie node. */
static igraph_error_t igraph_i_trie_init_node(igraph_trie_node_t *t) {
    IGRAPH_STRVECTOR_INIT_FINALLY(&t->strs, 0);
    IGRAPH_VECTOR_PTR_INIT_FINALLY(&t->children, 0);
    IGRAPH_VECTOR_INT_INIT_FINALLY(&t->values, 0);
    IGRAPH_FINALLY_CLEAN(3);
    return IGRAPH_SUCCESS;
}

static void igraph_i_trie_destroy_node(igraph_trie_node_t *t);

/**
 * \ingroup igraphtrie
 * \brief Creates a trie.
 *
 * \param t An uninitialized trie.
 * \param storekeys Specifies whether keys are stored for reverse lookup.
 * \return Error code: Errors by \ref igraph_strvector_init(),
 *         \ref igraph_vector_ptr_init() and \ref igraph_vector_init() might be returned.
 */

igraph_error_t igraph_trie_init(igraph_trie_t *t, igraph_bool_t storekeys) {
    t->maxvalue = -1;
    t->storekeys = storekeys;
    IGRAPH_CHECK(igraph_i_trie_init_node(&t->node));
    IGRAPH_FINALLY(igraph_i_trie_destroy_node, &t->node);
    if (storekeys) {
        IGRAPH_CHECK(igraph_strvector_init(&t->keys, 0));
    }

    IGRAPH_FINALLY_CLEAN(1);
    return IGRAPH_SUCCESS;
}

static void igraph_i_trie_destroy_node_helper(igraph_trie_node_t *t, igraph_bool_t sfree) {
    igraph_strvector_destroy(&t->strs);
    igraph_integer_t children_size = igraph_vector_ptr_size(&t->children);
    for (igraph_integer_t i = 0; i < children_size; i++) {
        igraph_trie_node_t *child = VECTOR(t->children)[i];
        if (child != NULL) {
            igraph_i_trie_destroy_node_helper(child, true);
        }
    }
    igraph_vector_ptr_destroy(&t->children);
    igraph_vector_int_destroy(&t->values);
    if (sfree) {
        IGRAPH_FREE(t);
    }
}

/* Deallocates a trie node. */
static void igraph_i_trie_destroy_node(igraph_trie_node_t *t) {
    igraph_i_trie_destroy_node_helper(t, false);
}

/**
 * \ingroup igraphtrie
 * \brief Destroys a trie (frees allocated memory).
 *
 * \param t The trie.
 */

void igraph_trie_destroy(igraph_trie_t *t) {
    if (t->storekeys) {
        igraph_strvector_destroy(&t->keys);
    }
    igraph_i_trie_destroy_node(&t->node);
}


/* Computes the location (index) of the first difference between 'str' and 'key' */
static size_t igraph_i_strdiff(const char *str, const char *key) {
    size_t diff = 0;
    while (key[diff] != '\0' && str[diff] != '\0' && str[diff] == key[diff]) {
        diff++;
    }
    return diff;
}

/**
 * \ingroup igraphtrie
 * \brief Search/insert in a trie (not to be called directly).
 *
 * \return Error code, usually \c IGRAPH_ENOMEM.
 */

static igraph_error_t igraph_i_trie_get_node(
    igraph_trie_node_t *t, const char *key, igraph_integer_t newvalue,
    igraph_integer_t *id
) {
    assert(key != NULL);

    /* If newvalue is negative, we don't add the node if nonexistent, only check
     * for its existence */
    igraph_bool_t add = (newvalue >= 0);

    igraph_integer_t strs_size = igraph_strvector_size(&t->strs);
    for (igraph_integer_t i = 0; i < strs_size; i++) {
        size_t diff;
        const char *str = igraph_strvector_get(&t->strs, i);
        diff = igraph_i_strdiff(str, key);

        if (diff == 0) {

            /* ------------------------------------ */
            /* No match, next */

        } else if (str[diff] == '\0' && key[diff] == '\0') {

            /* ------------------------------------ */
            /* They are exactly the same */
            if (VECTOR(t->values)[i] != -1) {
                *id = VECTOR(t->values)[i];
                return IGRAPH_SUCCESS;
            } else {
                VECTOR(t->values)[i] = newvalue;
                *id = newvalue;
                return IGRAPH_SUCCESS;
            }

        } else if (str[diff] == '\0') {

            /* ------------------------------------ */
            /* str is prefix of key, follow its link if there is one */
            igraph_trie_node_t *node = VECTOR(t->children)[i];
            if (node != NULL) {
                return igraph_i_trie_get_node(node, key + diff, newvalue, id);
            } else if (add) {
                igraph_trie_node_t *new_node = IGRAPH_CALLOC(1, igraph_trie_node_t);
                IGRAPH_CHECK_OOM(new_node, "Cannot add to trie.");
                IGRAPH_FINALLY(igraph_free, new_node);

                IGRAPH_STRVECTOR_INIT_FINALLY(&new_node->strs, 1);
                IGRAPH_VECTOR_PTR_INIT_FINALLY(&new_node->children, 1);
                IGRAPH_VECTOR_INT_INIT_FINALLY(&new_node->values, 1);
                IGRAPH_CHECK(igraph_strvector_set(&new_node->strs, 0, key + diff));
                IGRAPH_FINALLY_CLEAN(4);

                VECTOR(new_node->children)[0] = 0;
                VECTOR(new_node->values)[0] = newvalue;

                VECTOR(t->children)[i] = new_node;

                *id = newvalue;
                return IGRAPH_SUCCESS;
            } else {
                *id = -1;
                return IGRAPH_SUCCESS;
            }

        } else if (key[diff] == '\0' && add) {

            /* ------------------------------------ */
            /* key is prefix of str, the node has to be cut */
            char *str2;

            igraph_trie_node_t *node = IGRAPH_CALLOC(1, igraph_trie_node_t);
            IGRAPH_CHECK_OOM(node, "Cannot add to trie.");
            IGRAPH_FINALLY(igraph_free, node);

            IGRAPH_STRVECTOR_INIT_FINALLY(&node->strs, 1);
            IGRAPH_VECTOR_PTR_INIT_FINALLY(&node->children, 1);
            IGRAPH_VECTOR_INT_INIT_FINALLY(&node->values, 1);
            IGRAPH_CHECK(igraph_strvector_set(&node->strs, 0, str + diff));

            VECTOR(node->children)[0] = VECTOR(t->children)[i];
            VECTOR(node->values)[0] = VECTOR(t->values)[i];

            str2 = strdup(str);
            IGRAPH_CHECK_OOM(str2, "Cannot add to trie.");
            IGRAPH_FINALLY(igraph_free, str2);
            str2[diff] = '\0';

            IGRAPH_CHECK(igraph_strvector_set(&t->strs, i, str2));

            IGRAPH_FREE(str2);
            IGRAPH_FINALLY_CLEAN(5);

            VECTOR(t->values)[i] = newvalue;
            VECTOR(t->children)[i] = node;

            *id = newvalue;
            return IGRAPH_SUCCESS;

        } else if (add) {

            /* ------------------------------------ */
            /* the first diff characters match */
            char *str2;

            igraph_trie_node_t *node = IGRAPH_CALLOC(1, igraph_trie_node_t);
            IGRAPH_CHECK_OOM(node, "Cannot add to trie.");
            IGRAPH_FINALLY(igraph_free, node);

            IGRAPH_STRVECTOR_INIT_FINALLY(&node->strs, 2);
            IGRAPH_VECTOR_PTR_INIT_FINALLY(&node->children, 2);
            IGRAPH_VECTOR_INT_INIT_FINALLY(&node->values, 2);
            IGRAPH_CHECK(igraph_strvector_set(&node->strs, 0, str + diff));
            IGRAPH_CHECK(igraph_strvector_set(&node->strs, 1, key + diff));
            VECTOR(node->children)[0] = VECTOR(t->children)[i];
            VECTOR(node->children)[1] = 0;
            VECTOR(node->values)[0] = VECTOR(t->values)[i];
            VECTOR(node->values)[1] = newvalue;

            str2 = strdup(str);
            IGRAPH_CHECK_OOM(str2, "Cannot add to trie.");

            str2[diff] = '\0';
            IGRAPH_FINALLY(igraph_free, str2);

            IGRAPH_CHECK(igraph_strvector_set(&t->strs, i, str2));

            IGRAPH_FREE(str2);
            IGRAPH_FINALLY_CLEAN(5);

            VECTOR(t->values)[i] = -1;
            VECTOR(t->children)[i] = node;

            *id = newvalue;
            return IGRAPH_SUCCESS;
        } else {

            /* ------------------------------------------------- */
            /* No match, but we requested not to add the new key */
            *id = -1;
            return IGRAPH_SUCCESS;
        }
    }

    /* ------------------------------------ */
    /* Nothing matches */

    if (add) {
        /* Memory saving at the cost of performance may be possible by using the pattern
         *     CHECK(reserve(vec, size(vec) + 1));
         *     push_back(vec, value);
         * This was the original pattern used before igraph 0.10. */
        IGRAPH_CHECK(igraph_strvector_push_back(&t->strs, key));
        IGRAPH_CHECK(igraph_vector_ptr_push_back(&t->children, NULL));
        IGRAPH_CHECK(igraph_vector_int_push_back(&t->values, newvalue));
        *id = newvalue;
    } else {
        *id = -1;
    }

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup igraphtrie
 * \brief Search/insert a null-terminated string in a trie.
 *
 * \param t The trie.
 * \param key The string to search for. If not found, it will be inserted.
 * \param id The index of the string is stored here.
 * \return Error code, usually \c IGRAPH_ENOMEM.
 */

igraph_error_t igraph_trie_get(igraph_trie_t *t, const char *key, igraph_integer_t *id) {
    assert(key != NULL);

    if (*key == '\0') {
        IGRAPH_ERROR("Keys in a trie cannot be empty.", IGRAPH_EINVAL);
    }

    if (!t->storekeys) {
        IGRAPH_CHECK(igraph_i_trie_get_node(&t->node, key, t->maxvalue + 1, id));
        if (*id > t->maxvalue) {
            t->maxvalue = *id;
        }
    } else {
        igraph_error_t ret;

        IGRAPH_FINALLY_ENTER();
        /* Add it to the string vector first, we can undo this later */
        ret = igraph_strvector_push_back(&t->keys, key);
        if (ret != IGRAPH_SUCCESS) {
            IGRAPH_FINALLY_EXIT();
            IGRAPH_ERROR("Cannot get element from trie.", ret);
        }
        ret = igraph_i_trie_get_node(&t->node, key, t->maxvalue + 1, id);
        if (ret != IGRAPH_SUCCESS) {
            igraph_strvector_resize(&t->keys, igraph_strvector_size(&t->keys) - 1); /* shrinks, error safe */
            IGRAPH_FINALLY_EXIT();
            IGRAPH_ERROR("Cannot get element from trie.", ret);
        }

        /* everything is fine */
        if (*id > t->maxvalue) {
            t->maxvalue = *id;
        } else {
            igraph_strvector_resize(&t->keys, igraph_strvector_size(&t->keys) - 1); /* shrinks, error safe */
        }
        IGRAPH_FINALLY_EXIT();
    }

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup igraphtrie
 * \brief Search/insert a string of given length in a trie.
 *
 * This function is identical to \ref igraph_trie_get(), except that
 * it takes a string of a given length as input instead of a null-terminated
 * string.
 *
 * \param t The trie.
 * \param key The string to search for. If not found, it will be inserted.
 * \param length The length of \p key.
 * \param id The index of the string is stored here.
 * \return Error code, usually \c IGRAPH_ENOMEM.
 */

igraph_error_t igraph_trie_get_len(
        igraph_trie_t *t, const char *key,
        igraph_integer_t length,
        igraph_integer_t *id) {

    char *tmp = strndup(key, length);
    IGRAPH_CHECK_OOM(tmp, "Cannot get from trie.");
    IGRAPH_FINALLY(igraph_free, tmp);
    IGRAPH_CHECK(igraph_trie_get(t, tmp, id));
    IGRAPH_FREE(tmp);
    IGRAPH_FINALLY_CLEAN(1);

    return IGRAPH_SUCCESS;
}

/**
 * \ingroup igraphtrie
 * \brief Search in a trie.
 *
 * This variant does not add \p key to the trie if it does not exist.
 * In this case, a negative \p id is returned.
 *
 * \param t The trie.
 * \param key The string to search for.
 * \param id If \p key is found, its index is stored here. Otherwise,
 *    a negative value is returned.
 * \param Error code.
 */

igraph_error_t igraph_trie_check(igraph_trie_t *t, const char *key, igraph_integer_t *id) {
    IGRAPH_CHECK(igraph_i_trie_get_node(&t->node, key, -1, id));
    return IGRAPH_SUCCESS;
}

/**
 * \ingroup igraphtrie
 * \brief Get an element of a trie based on its index.
 *
 * \param t The trie.
 * \param idx The index of the string. It is not checked that it is within range.
 * \return The string with the given index. If the trie does not store the keys for
 *   reverse lookup, \c NULL is returned.
 */

const char* igraph_trie_idx(igraph_trie_t *t, igraph_integer_t idx) {
    if (! t->storekeys) {
        return NULL;
    }
    return igraph_strvector_get(&t->keys, idx);
}

/**
 * \ingroup igraphtrie
 * \brief Returns the size of a trie.
 *
 * \param t The trie.
 * \return The size of the trie, i.e. one larger than the maximum index.
 */

igraph_integer_t igraph_trie_size(igraph_trie_t *t) {
    return t->maxvalue + 1;
}

/* Hmmm, very dirty.... */

/**
 * \ingroup igraphtrie
 * \brief Retrieves all the keys from the trie.
 *
 * </para><para>
 * Note that the returned pointer is a \em borrowed reference into the internal
 * string vector of the trie. Do \em not modify it and do \em not use it after
 * the trie was destroyed.
 *
 * \param t The trie.
 * \return The borrowed reference.
 */

const igraph_strvector_t* igraph_i_trie_borrow_keys(igraph_trie_t *t) {
    return &t->keys;
}
