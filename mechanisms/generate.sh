for mech in pas hh
do
    ../external/modparser/bin/modcc -t cpu -o ../include/mechanisms/$mech.hpp ./mod/$mech.mod
done
