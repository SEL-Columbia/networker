from networker.classes.unionfind import PriorityQueue
from nose.tools import eq_


def test_push_pop():
    """
    Tests that pushing a list of items
    into a queue and then popping them out
    results in the expected behavior.
    """

    q = PriorityQueue()

    # input list (obj, priority) should be reversed
    # in the priority_queue
    input_list = [((1), 9), ((2), 8), ((3), 7),
                  ((4), 6), ((5), 5), ((6), 4),
                  ((7), 3), ((8), 2), ((9), 1)]

    # insert the items in the queue
    for obj, p in input_list:
        q.push(obj, p)

    # pop the items into another list
    output = []
    while q._queue:
        output.append(q.pop())

    # make sure it lines up with expected result
    eq_(output, range(1, 10)[::-1])


def test_peek_top():
    """
    Tests that peeking at the top of the priortity queue
    leaves it unmodified
    """

    q = PriorityQueue()

    # input list (obj, priority) should be reversed
    # in the priority_queue
    input_list = [((1), 9), ((2), 8), ((3), 7),
                  ((4), 6), ((5), 5), ((6), 4),
                  ((7), 3), ((8), 2), ((9), 1)]

    # insert the items in the queue
    for obj, p in input_list:
        q.push(obj, p)

    # save the state of the queue
    state = [i for i in q._queue]

    # peek a few times
    [q.top() for i in range(100)]

    eq_(state, q._queue)


def test_merge():
    """
    Tests that merging two priority queues
    results in the expected behavior
    """

    q1 = PriorityQueue()
    q2 = PriorityQueue()

    # input list (obj, priority) should be reversed
    # in the priority_queue
    input_list = [((1), 9), ((2), 8), ((3), 7),
                  ((4), 6), ((5), 5), ((6), 4),
                  ((7), 3), ((8), 2), ((9), 1)]

    # insert the items into the queues
    for idx, (obj, p) in enumerate(input_list):
        if idx < 4:
            q1.push(obj, p)
        else:
            q2.push(obj, p)

    # merge the queues
    q1.merge(q2)

    # dump the queue into a list
    output = []
    while q1._queue:
        output.append(q1.pop())

    # validate the expected behavior
    eq_(output, range(1, 10)[::-1])
