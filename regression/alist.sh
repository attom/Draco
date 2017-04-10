#!/bin/bash

#------------------------------------------------------------------------------#
# Generate a list of authors by using git annotate
#------------------------------------------------------------------------------#

# Should be run from the draco top-level source directory.
if ! test -f CMakeLists.txt; then
  echo "ERROR: Must be run from the top level draco source directory."
  exit 1
fi
if ! test -f README.md; then
  echo "ERROR: Must be run from the top level draco source directory."
  exit 1
fi

run () {
   echo $1
   if ! test $dry_run; then
      eval $1
   fi
}

function math ()
{
    calc="${@//p/+}";
    calc="${calc//x/*}";
    bc -l <<< "scale=10;$calc"
}

#------------------------------------------------------------------------------#
# List of contributors
author_list_file=`mktemp`
git shortlog -s | cut -c8- &> $author_list_file

# Put the user list into an array (some entries may have spaces).
IFS=$'\r\n' GLOBIGNORE='*' command eval 'user_list=($(cat $author_list_file))'

#------------------------------------------------------------------------------#
# http://stackoverflow.com/questions/1265040/how-to-count-total-lines-changed-by-a-specific-author-in-a-git-repository

# Ingore some files used in Jayenne:
# ignorespec=":(exclude)*/output.in" ":(exclude)*/*.bench"

# count loc for each author and report.
author_loc=`mktemp`
for name in "${user_list[@]}"; do

  if [[ $name == "regress" ]]; then
    continue
  fi

  # sort by net lines added (lines added - lines removed)
  numlines=`git log --author="$name" --pretty=tformat: --numstat -- . ":(exclude)clubimc" ":(exclude)wedgehog" ":(exclude)milagro" ":(exclude)*/output.in" ":(exclude)*/*.bench" | awk '{ add += $1; subs += $2; loc += $1 - $2; sum += $1 + $2 } END { printf "%s:\n", loc}'`

  merge_name=$name
  case $name in
    along | Alex*  ) merge_name="Alex Long" ;;
    clevelam) merge_name="Matt Cleveland" ;;
    gaber) merge_name="Gabe Rockefeller" ;;
    jdd) merge_name="Jeff Densmore" ;;
    jhchang | Jae* ) merge_name="Jae Chang" ;;
    keadyk | Kendra* ) merge_name="Kendra P. Keady" ;;
    kellyt | kgt | 107638 | Kelly*) merge_name="Kelly Thompson" ;;
    kgbudge) merge_name="Kent G. Budge" ;;
    kwang) merge_name="Katherine Wang" ;;
    lowrie) merge_name="Rob Lowrie" ;;
    lpritch) merge_name="Lori Pritchett-Sheats" ;;
    maxrosa) merge_name="Massimiliano Rosa" ;;
    ntmyers) merge_name="Nick Meyers" ;;
    pahrens) merge_name="Peter Ahrens" ;;
    talbotp) merge_name="Paul Talbot" ;;
    tmonster) merge_name="Todd Urbatsch" ;;
    warsa) merge_name="James Warsa" ;;
    wollaber) merge_name="Allan Wollaber" ;;
    wollaege) merge_name="Ryan Thomas Wollaeger" ;;
  esac

  #echo "$numlines$merge_name"
  echo "$merge_name:$numlines" | sed -e 's/:$//'

done | sort -rn > $author_loc

# cat $author_loc

# Merge data

prev_author="none"
prev_loc=0

# Put the user list into an array (some entries may have spaces).
IFS=$'\r\n' GLOBIGNORE='*' command eval 'entries=($(cat $author_loc))'
rm $author_loc
for entry in "${entries[@]}"; do

  current_author=`echo $entry | sed -e 's/:.*//'`
  current_loc=`echo $entry | sed -e 's/.*://'`
  if [[ $current_loc == "" ]]; then current_loc=0; fi
  if [[ $current_author == $prev_author ]]; then
    prev_loc=`math "$prev_loc + $current_loc"`
  else
    # previous author information before starting a new author.
    if ! [[ $prev_author == "none" ]]; then
      echo "${prev_loc}:${prev_author}" >> $author_loc
    fi

    prev_author=$current_author
    prev_loc=$current_loc
  fi

done
# the last entry
echo "${prev_loc}:${prev_author}" >> $author_loc

cat $author_loc | sort -rn

# cleanup
rm $author_loc $author_list_file

# while read line; do
#   echo $line | sed -e 's/\([0-9]*\):/current_developers[\1]=/'
# done < $author_loc

#------------------------------------------------------------------------------#
# Old svn based script modified for git
# # Build a list of files to inspect.
# files=`find . -name '*.hh' -o -name '*.cc' -o -name '*.txt' \
#          -o -name '*.cmake' -o -name '*.in' -o -name '*.h'`

# # generate full annotations for these files.
# author_per_loc=`mktemp`
# for file in $files; do
#   echo "annotate $file"
#   git annotate --line-porcelain $file |  grep '^author ' | sed -e 's/author //' >> $author_per_loc
# done

# # generate a list of users identified by 'git annotate'
# #user_list=`cat $author_per_loc | awk '{print $3}' | sort -u`
# author_list=`mktemp`
# cat $author_per_loc | sort -u | grep -v Not > $author_list

# # Put the user list into an array (some entries may have spaces).
# IFS=$'\r\n' GLOBIGNORE='*' command eval 'user_list=($(cat $author_list))'

# # count loc for each author and report.
# author_loc=`mktemp`
# for name in "${user_list[@]}"; do
#   if ! test "$name" == "kgt" && ! test "$name" == "Kelly (KT) Thompson" && \
#      ! test "$name" == 107638; then
#     numlines=`grep -c "$name" $author_per_loc`
#     echo "$numlines:$name"
#   fi
# done | sort -rn > $author_loc

# while read line; do
#   echo $line | sed -e 's/\([0-9]*\):/current_developers[\1]=/'
# done < $author_loc



#------------------------------------------------------------------------------#
# end alist.sh
#------------------------------------------------------------------------------#
