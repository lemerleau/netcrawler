#!/bin/bash
tr="1";
for id in `ps -aux | grep 'user_id'| cut -f 2 -d ' '` ; 

do 
	kill -9 $id
done
