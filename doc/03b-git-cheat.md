\cond HIDDEN

# git tips

* Mostly for

## git diff

```shell
git diff --stat old_branch new_branch
```

git diff --stat old_branch
 (same as git diff --stat old_branch current_branch)

## git bisect

```bash
git bisect start
git bisect bad dev
git bisect good b36e8dba

git bisect run bash -lc '
  rm Makefile
  cp ./doc/Makefile .              2>/dev/null || :
  cp ./doc/examples/Makefile .     2>/dev/null || :
  make clean
  make -j12 tests
  ./tests <specific test>
'
```

\endcond
