require File.expand_path(__FILE__).sub(%r(/test/.*), '/test/test_helper.rb')
#require File.expand_path(__FILE__).sub(%r(.*/test/), '').sub(/test_(.*)\.rb/,'\1')

require 'scout'
class TestClass < Test::Unit::TestCase
  def self.collapse_offsets(offsets)
    last = nil
    start = nil
    segments = []
    current = []
    while i = offsets.shift
      if current.any? && current.last + 1 == i
        current << i
      elsif current.any?
        segments << [current.first, current.last]
        current = []
      else
        current = []
      end
      current << i
    end
    segments << [current.first, current.last] if current.any?

    puts segments.inspect
    #segments.collect{|s,e| [s,e].uniq.collect{|i| AGS::TIME_POINTS[i] } * "-" }
  end

  def test_collapse
    TestClass.collapse_offsets([1,4])
  end
end

