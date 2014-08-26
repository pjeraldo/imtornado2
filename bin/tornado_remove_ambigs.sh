#!/usr/bin/env bash
#Flatten, then match for an ambig, remove said line, and the corresponding id line
#Much faster than python for this task
perl -pe 's/$/xXx/ if /^>/; s/^/xXx/ if /^>/' $1 |tr -d '\n'|perl -pe 's/xXx/\n/g' |sed '/^$/d'| sed -n '/[WSMKRYBDHVNwsmkrybdhvn]/{/>/!{s/.*//;x;d;};};x;p;${x;p;}'| sed '/^$/d' > $2
