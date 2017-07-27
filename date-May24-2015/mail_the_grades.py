#!/usr/bin/python

import sys
import os
import re
import time
import smtplib
from email.mime.text import MIMEText

#-----------------------------------------------------------------------------
# Student
#-----------------------------------------------------------------------------
class Student:
    """Student's info and complete grade"""
    def __init__(self, email) :
        self.email = email
    name = ''
    rank = ''
    email = ''
    grade = 0
    breakdown = ''
    comments = ''

    def assign_grade(self, breakdown) :
        if len(breakdown) == 0:
            return
        _pattern = re.compile('\s|(?<!\d)[,.]|[,.](?!\d)')
        gparts = _pattern.split(breakdown)
        if len(gparts) == 0 :
            return
        for g in gparts :
            if len(g) == 0:
                continue
            part = int(g)
            self.grade += part
            self.breakdown += '%d, ' % part
        self.breakdown = self.breakdown.rstrip(', ') 

    def __str__(self) :
        s = "Name: %s (%s)\n" \
        "e-mail: %s\n" \
        "Grade total: %d\n" \
        "Detailed breakdown: %s\n" \
        "Grading comments: %s\n" % (self.name, self.rank, self.email,
                                    self.grade, self.breakdown,
                                    self.comments)
        s += '-----------------------------------\n'
        return s

#-----------------------------------------------------------------------------
# content
#-----------------------------------------------------------------------------
def prep_email(csc, assignment, intro, student) :
    content  = 'Dear %s,\n' % student.name
    content += '\n'    
    content += 'Below is your grade for %s\n' % assignment
    content += 'and some comments regarding the grading for this assignment.\n'
    content += '(Comments are optional --- you might not get any, especially \n'
    content += 'for sections you got full grade.) \n'
    content += '\n'        
    content += intro
    content += '\n'
    content += 'Your grade for this assignment is:  %d \n' % student.grade
    content += 'broken down as follows, according to the previous:\n%s\n' % \
               student.breakdown
    content += '\n'
    content += 'A few comments regarding your grade follow below:'
    content += '\n\n'
    content += student.comments
    content += '\n\n'

    subject = '%s --- %s' % (csc, assignment)

    return (subject, content)
    
#-----------------------------------------------------------------------------
# email
#-----------------------------------------------------------------------------
def send_email(from_addr, to_addr, subject, content):

    msg = MIMEText(content)
    msg['Subject'] = subject
    msg['From'] = from_addr
    msg['To'] = to_addr

    print ''
    print '-------------------------------------------------------------'
    print ''

    # print subject
    # print from_addr
    # print to_addr
    # print content

    # print ''
    # print '-------------------------------------------------------------'
    # print ''

    print '*** Sending e-mail to %s *** ' % to_addr
    time.sleep(1)
    s = smtplib.SMTP('localhost')
    s.sendmail(from_addr, [to_addr], msg.as_string())

    return

#-----------------------------------------------------------------------------
# process_file
#-----------------------------------------------------------------------------
def process_file(fname) :
    
    lines = None

    try :
        fd = open(fname, 'r')
    except IOError :
        print 'ERROR: Cannot open file %s'  % fname
        sys.exit(1)
    else :
        lines = fd.readlines()
        fd.close()

    csc = ''
    from_addr = ''
    assignment = ''
    intro = ''
    students = []


    for l in lines :
        if l.startswith('#') :
            appendto = ''
        if appendto == 'intro' :
            intro += l
            continue
        if appendto == 'comments' :
            students[0].comments += l
            continue
        if l.startswith('#-') :
            continue
        if l.startswith('#ClassName:') :
            param, csc = l.split(':')
            csc = csc.strip()
            continue
        if l.startswith('#ClassEmail:') :
            param, from_addr = l.split(':')
            from_addr = from_addr.strip()
            continue
        if l.startswith('#Assignment:') :
            param, assignment = l.split(':')
            assignment = assignment.strip()
            continue
        if l.startswith('#GradingCommon:') :
            param, intro = l.split(':')
            intro = intro.strip()
            appendto = 'intro'
            continue
        if l.startswith('#StudentEmail:') :
            param, email = l.split(':')
            email = email.strip()
            student = Student(email)
            students.insert(0, student)
            continue
        if l.startswith('#StudentName:') :
            param, name = l.split(':')
            name = name.strip()
            students[0].name = name
            continue
        if l.startswith('#Rank:') :
            param, rank = l.split(':')
            rank = rank.strip()
            students[0].rank = rank
            continue
        if l.startswith('#GradingBreakdown:') :
            param, breakdown = l.split(':')
            breakdown = breakdown.strip()
            students[0].assign_grade(breakdown)
            continue
        if l.startswith('#GradingComments:') :
            param, comments = l.split(':')
            comments = comments.strip()
            students[0].comments = comments
            appendto = 'comments'            
            continue

    return from_addr, csc, assignment, intro, reversed(students)

    return

#-----------------------------------------------------------------------------
# main
#-----------------------------------------------------------------------------
def main(args):

    if len(sys.argv) != 2 :
        print 'USAGE: %s grades-file.txt ' % args[0]
        print 'where:'
        print '       grades-file.txt is a properly formatted file;'
        print '       see grades-template.txt or previous assignments'
        print '       for an example'
        sys.exit(1)

    (from_addr, csc, assignment, intro, students) = process_file(args[1])

    for student in students :
        subject, content = prep_email(csc, assignment, intro, student)
        send_email(from_addr, student.email, subject, content)

    return

#-----------------------------------------------------------------------------
# script entry point
#-----------------------------------------------------------------------------
if __name__ == '__main__' :
    main(sys.argv)
