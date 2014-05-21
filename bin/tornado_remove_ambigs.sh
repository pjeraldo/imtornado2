#!/usr/bin/env bash
#Flatten, then match for an ambig, remove said line, and the corresponding id line
#Much faster than python for this task
sed '/>/{s/^/xXx/;s/$/xXx/;}' $1 |tr -d '\n'|sed 's/xXx/\n/g' |sed '/^$/d'| sed -n '/[WSMKRYBDHVNwsmkrybdhvn]/{/>/!{s/.*//;x;d;};};x;p;${x;p;}'| sed '/^$/d' > $2
