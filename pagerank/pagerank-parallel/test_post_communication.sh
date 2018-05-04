sort test_outgoing_to_incoming_main.txt > test1.txt
sort test_outgoing_to_incoming_post_comm.txt > test2.txt
diff test1.txt test2.txt
