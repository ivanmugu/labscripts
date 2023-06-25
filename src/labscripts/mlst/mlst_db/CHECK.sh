#!/bin/sh
#
# Checks the consistency of the MLST database.
# Errors are written to stdout.

err_msg() {
    echo "[ERROR] $*"
}

warn_msg() {
    echo "[WARN]  $*"
}

# Check that every scheme in config has FSA and TSV
grep -E '^[^#]' config | while read SCHEME REST; do
    FSA_FILE="${SCHEME}/${SCHEME}.fsa"
    TSV_FILE="${SCHEME}/${SCHEME}.tsv"
    [ -f "$FSA_FILE" ] || err_msg "not found: $FSA_FILE"
    [ -f "$TSV_FILE" ] || err_msg "not found: $TSV_FILE"
done

# Check that every FSA has a TSV
find . -name '*.fsa' | while read FSA_FILE; do
    TSV_FILE="${FSA_FILE%.fsa}.tsv"
    [ -f "$TSV_FILE" ] || err_msg "no TSV for: $FSA_FILE"
done

# Check that every TSV has an FSA
find . -name '*.tsv' | while read TSV_FILE; do
    FSA_FILE="${TSV_FILE%.tsv}.fsa"
    [ -f "$FSA_FILE" ] || err_msg "no FSA for: $TSV_FILE"
done

# Check that every FSA is in config
find . -name '*.fsa' | while read FSA_FILE; do
    SCHEME="$(basename "$FSA_FILE" ".fsa")"
    cut -f1 config | grep -q "^${SCHEME}\$" || err_msg "not in config: $SCHEME"
done

# Check that the genes listed in config are in the TSV headers (in any order)
grep -E '^[^#].+$' config | cut -f1,3 | while read SCHEME GENES; do
    unset MISSING
    echo "$GENES" | tr ',' ' ' | while read GENE; do
        head -1 "$SCHEME/$SCHEME.tsv" | grep -xqE "(.*\t)?${GENE}(\t.*)?" || MISSING="$MISSING $GENE"
    done
    [ -z "$MISSING" ] || err_msg "missing gene(s) in TSV file for $SCHEME: $MISSING"
    # (Only) if none were missing, also check that order is the same as in config
    [ -n "$MISSING" ] || head -1 "$SCHEME/$SCHEME.tsv" | tr '\t' ',' | grep -xqE "(.*,)?${GENES}(,.*)?" ||
        warn_msg "order of genes in TSV different from config for $SCHEME"
done

# Check that the genes in the FSA match those in the config
grep -E '^[^#].+$' config | cut -f1,3 | while read SCHEME GENES; do
    FSA_GENES="$(grep '^>' "$SCHEME/$SCHEME.fsa" | sed -Ee 's/>(.*)[^0-9][0-9]+( .*)?$/\1/' | sort -u | tr '\n' ',')"
    CFG_GENES="$(echo $GENES | tr ',' '\n' | sort | tr '\n' ',')"
    [ "$FSA_GENES" = "$CFG_GENES" ] || err_msg "genes in FSA don't match config for $SCHEME: $FSA_GENES vs $CFG_GENES"
done

# vim: sts=4:sw=4:si:ai:et
